"""
This makes a digital surface model (DSM), digital elevation model (DEM), and
canopy height model (CHM) from a LAZ file point cloud produced by structure from
motion processing of overlapping drone imagery by Pix4DMapper.
"""

import os
import sys
import numpy
import argparse
from numba import jit
from osgeo import gdal
from osgeo import osr
from scipy import ndimage
from pylidar import lidarprocessor
from pylidar.lidarformats import generic
from pylidar.toolbox import interpolation
from pylidar.toolbox.grdfilters import pmf


def process_laz(lazFile, epsg, pixelSize, projName, outDir, maxWinSize):
    """
    Main function to create a digital surface model, digital elevation model,
    and canopy height model.
    """
    
    # Create names for output files
    dsmFile = os.path.join(outDir, r'%s_dsm.tif'%projName)
    demFile = os.path.join(outDir, r'%s_dem.tif'%projName)
    chmFile = os.path.join(outDir, r'%s_chm.tif'%projName)
    
    # Read in the point cloud and calculate image extents
    points = readLidarPoints(lazFile, colNames=['X', 'Y', 'Z'])
    (x, y, z) = (points['X'], points['Y'], points['Z'])
    minX, maxX, minY, maxY, minZ, maxZ = get_mmXYZ(x,y,z)
    if minX == int(minX):
        minX = minX - 1
    if minY == int(minY):
        minY = minY - 1
    tileXsize = numpy.ceil(maxX) - int(minX)
    tileYsize = numpy.ceil(maxY) - int(minY)
    nRows = int(numpy.ceil(tileYsize / pixelSize))
    nCols = int(numpy.ceil(tileXsize / pixelSize))   
    info = generic.getLidarFileInfo(lazFile)
    inFormat = info.getDriverName()
    (row, col) = xyToRowCol(x, y, int(minX), int(numpy.ceil(maxY)), pixelSize)

    # Create digital surface model
    nullVal = -999
    maxHeight = numpy.zeros((nRows, nCols), dtype=numpy.float32) + nullVal
    maxGridding(maxHeight, row, col, z)
    nullArray = (maxHeight == nullVal)

    pxlCoords = get_uneven_grid(nCols, nRows, int(minX), int(numpy.ceil(maxY)),
                                pixelSize)
    xVals = pxlCoords[0][~nullArray]
    yVals = pxlCoords[1][~nullArray]
    zVals = maxHeight[~nullArray]
    dsm = interpolation.interpGrid(xVals, yVals, zVals, pxlCoords, 'pynn')
    dsm[numpy.isnan(dsm)] = nullVal
    dsm[~nullArray] = maxHeight[~nullArray]
    dsm = numpy.round(dsm.astype(numpy.float32), 3)
    writeImage(dsm, dsmFile, driver='GTiff', tlx=int(minX),
               tly=int(numpy.ceil(maxY)), binsize=pixelSize, epsg=epsg,
               nullVal=nullVal)

    # Create digital elevation model
    nullVal = 9999.0
    minHeight = numpy.zeros((nRows, nCols), dtype=numpy.float32) + nullVal
    minGridding(minHeight, row, col, z)
    nullArray = (minHeight == nullVal)
     
    ############################################################################
    # This section implements a progressive morphological filter, modified from
    # the pmf.applyPMF() function in pylidar.
    dataArr = minHeight
    noDataMask = ~nullArray
    binGeoSize = pixelSize
    initWinSize = 1
    
    # Keep number of filters to ~10 by changing winSizeInc
    winSizeInc = int(maxWinSize/10.0)

    slope = 0.3
    dh0 = 0.3
    dhmax = 5
    k = numpy.arange(0, maxWinSize, winSizeInc)
    winSize = (2*k*initWinSize) + 1
    A = dataArr
    A = A.astype(numpy.float64)
    nCols = A.shape[0]
    nRows = A.shape[1]
    A = pmf.doNearestNeighbourInterp(A, noDataMask, nCols, nRows)
    winSize_tminus1 = numpy.zeros([winSize.shape[0]])
    winSize_tminus1[1:] = winSize[:-1]
    A = pmf.doOpening(A, maxWinSize, winSize_tminus1, binGeoSize, slope, dh0,
                      dhmax)
    dem = ndimage.morphology.grey_dilation(A.astype(numpy.float32), size=(3,3))
    ############################################################################
    
    nullArray = (dem == nullVal)
    xVals = pxlCoords[0][~nullArray]
    yVals = pxlCoords[1][~nullArray]
    zVals = dem[~nullArray]
    dem = interpolation.interpGrid(xVals, yVals, zVals, pxlCoords, 'pynn')
    invalid = numpy.isnan(dem)
    nullVal = -999
    dem[invalid] = nullVal
    dem = numpy.round(dem.astype(numpy.float32), 3)
    writeImage(dem, demFile, driver='GTiff', tlx=int(minX),
               tly=int(numpy.ceil(maxY)), binsize=pixelSize, epsg=epsg,
               nullVal=nullVal)

    # Create canopy height model
    chm = dsm - dem
    chm[chm < 0] = 0
    chm[invalid] = nullVal
    writeImage(chm, chmFile, driver='GTiff', tlx=int(minX),
               tly=int(numpy.ceil(maxY)), binsize=pixelSize, epsg=epsg,
               nullVal=nullVal)
    
    print("Processing completed")


@jit
def maxGridding(grid, row, col, prop):
    """
    Create grid of maximum value of prop.
    """
    numPts = len(row)
    (nRows, nCols) = grid.shape
    for i in range(numPts):
        (r, c) = (row[i], col[i])
        if prop[i] > grid[r, c]:
            grid[r, c] = prop[i]


@jit
def minGridding(grid, row, col, prop):
    """
    Create grid of minimum value of prop.
    """
    numPts = len(row)
    (nRows, nCols) = grid.shape
    for i in range(numPts):
        (r, c) = (row[i], col[i])
        if prop[i] < grid[r, c]:
            grid[r, c] = prop[i]                             


def disk(radius):
    """
    Generates a flat, disk-shaped structuring element.
    """
    L = numpy.arange(-radius, radius + 1)
    X, Y = numpy.meshgrid(L, L)
    return numpy.array((X ** 2 + Y ** 2) <= radius ** 2, dtype=numpy.uint8)

            
def writeImage(image, outfile, driver='GTiff', tlx=0.0, tly=0.0, binsize=0.0,
               epsg=None, nullVal=None):
    """
    Write data to a GDAL supported image file format
    """
    if len(image.shape)==2:
        ny,nx = image.shape
        nz=1
    if len(image.shape)==3:
        nz,ny,nx = image.shape            
    driver = gdal.GetDriverByName(driver)
    dt = image.dtype
    
    if dt == 'uint8': gdaldtype = gdal.GDT_Byte
    if dt == 'int16': gdaldtype = gdal.GDT_Int16
    if dt == 'uint16': gdaldtype = gdal.GDT_UInt16
    if dt == 'int32': gdaldtype = gdal.GDT_Int32
    if dt == 'float32': gdaldtype = gdal.GDT_Float32
    if dt == 'float64': gdaldtype = gdal.GDT_Float64
    
    ds = driver.Create(outfile, nx, ny, nz, gdaldtype, ['COMPRESS=LZW'])
    ds.SetGeoTransform([tlx,binsize,0,tly,0,-binsize])

    if epsg is not None:
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(epsg)
        ds.SetProjection(proj.ExportToWkt())
    if nz>1:
        for i in range(nz):
            band = ds.GetRasterBand(i+1)
            band.WriteArray(image[i,:,:],0,0)
    else:
        band = ds.GetRasterBand(1)    
        band.WriteArray(image,0,0)
    
    # Set the null value on every band
    if nullVal is not None:
        for i in range(nz):
            band = ds.GetRasterBand(i+1)
            band.SetNoDataValue(nullVal)

    ds.FlushCache()
    

def get_mmXYZ(x,y,z):
    """
    basic XYZ min/max extraction
    """
    minX=numpy.min(x)
    maxX=numpy.max(x)
    minY=numpy.min(y)
    maxY=numpy.max(y)
    minZ=numpy.min(z)
    maxZ=numpy.max(z)
    return minX,maxX,minY,maxY,minZ,maxZ
    

def get_uneven_grid(xr, yr, xst, yst, psize):
    """
    Builds paired x and y grid locations for interpolation with pynninterp
    """  
    x_id = (numpy.array(range(int(xr)))) * psize + xst
    x_ids = (0.5 * psize) + numpy.ones((int(yr), 1), numpy.float64) * x_id
    y_id = (numpy.floor(yst) - numpy.array(range(int(yr))) * psize)    
    y_ids = numpy.ones((int(xr), 1), numpy.float64) * y_id - (0.5 * psize)
    y_ids = numpy.rot90(y_ids, 3)
    pxlCoords = [x_ids, y_ids]
    return pxlCoords


def readLidarPoints(filename, groundOnly=False, boundingbox=None,
                    colNames=['X', 'Y', 'Z']):
    """
    Read the requested columns for the points in the given file, in a memory-
    efficient manner.
    """
    datafiles = lidarprocessor.DataFiles()
    datafiles.infile = lidarprocessor.LidarFile(filename, lidarprocessor.READ)
    otherargs = lidarprocessor.OtherArgs()
    otherargs.groundOnly = groundOnly
    otherargs.colNames = colNames
    otherargs.dataArrList = []
    otherargs.boundingbox = boundingbox
    controls = lidarprocessor.Controls()
    controls.setSpatialProcessing(False)
    lidarprocessor.doProcessing(selectColumns, datafiles, otherArgs=otherargs, controls=controls)
    nPts = sum([len(a) for a in otherargs.dataArrList])
    if nPts > 0:
        fullArr = numpy.zeros(nPts, dtype=otherargs.dataArrList[0].dtype)
        i = 0
        for dataArr in otherargs.dataArrList:
            numPts = len(dataArr)
            fullArr[i:i+numPts] = dataArr
            i += len(dataArr)
    else:
        fullArr = numpy.array([])
    return fullArr


def selectColumns(data, otherargs):
    """
    Called from pylidar's doProcessing. Read the next block of lidar points,
    select out the requested columns. 
    """
    dataArr = data.infile.getPoints(colNames=otherargs.colNames)

    if otherargs.groundOnly:
        if 'CLASSIFICATION' in otherargs.colNames:
            pntClass = data['CLASSIFICATION']
        else:
            pntClass = data.infile.getPoints(colNames='CLASSIFICATION')
        mask = (pntClass == generic.CLASSIFICATION_GROUND)
        dataArr = dataArr[mask]
    
    if otherargs.boundingbox is not None:
        (xmin, xmax, ymin, ymax) = otherargs.boundingbox
        (x, y) = (dataArr['X'], dataArr['Y'])
        mask = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))
        dataArr = dataArr[mask]

    if len(dataArr) > 0:
        # Stash in the list of arrays. 
        otherargs.dataArrList.append(dataArr)

        
def xyToRowCol(x, y, xMin, yMax, pixSize):
    """
    For the given pixel size and xMin, yMax, convert the given arrays of x and y
    into arrays of row and column in a regular grid across the tile.
    """
    col = ((x - numpy.floor(xMin)) / pixSize).astype(numpy.uint32)
    row = ((numpy.ceil(yMax) - y) / pixSize).astype(numpy.uint32)
    return (row, col)


def getCmdargs():
    """
    Get the command line arguments.
    """
    p = argparse.ArgumentParser(
            description=("This makes a digital surface model (DSM), digital "+
                         "elevation model (DEM), and canopy height model (CHM)"+
                         " from a LAZ file point cloud produced by structure "+
                         "from motion processing of overlapping drone imagery "+
                         "by Pix4DMapper."))
    p.add_argument("-i", "--inLaz", dest="inLaz", default=None,
                   help=("Input LAZ file"))
    p.add_argument("-e", "--epsg", dest="epsg", default=None,
                   help=("EPSG code for gridded outputs (e.g. 32755)"))
    p.add_argument("-p", "--projectName", dest="projectName", default=None,
                   help=("Project name for output file prefix"))
    p.add_argument("-o", "--outDir", dest="outDir", default=None,
                   help=("Directory for output images"))  
    p.add_argument("-s", "--pixelsize", dest="pixelsize", default=0.05,
                   help=("Pixelsize (m) for gridded outputs (default=0.05)"))
    p.add_argument("-w", "--windowsize", dest="windowsize", default=40,
                   help=("Maximum filter window size in pixels (default=40"))
    cmdargs = p.parse_args()
    if (cmdargs.inLaz is None or cmdargs.epsg is None or
        cmdargs.projectName is None or cmdargs.outDir is None):
        p.print_help()
        print("Must name input LAZ file, EPSG code, projectName and outDir.")
        sys.exit()
    return cmdargs


if __name__ == "__main__":
    cmdargs = getCmdargs()
    process_laz(cmdargs.inLaz, int(cmdargs.epsg), float(cmdargs.pixelsize),
                cmdargs.projectName, cmdargs.outDir, int(cmdargs.windowsize))
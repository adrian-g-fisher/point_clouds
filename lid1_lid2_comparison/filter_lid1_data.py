#!/usr/bin/env python

import os
import sys
import glob
import glob
import numpy
import laspy
from pyproj import crs
from alsveg import rw_als_data
from alsveg import gridding_interp_fns


def removeBrokenPulses(data):
    """
    Removes broken pulses from the input data recArray. Checks if
    TIMESTAMPs are the same and return_numbers are complete for each
    return in each pulse
    """
    firstReturnIndex = numpy.where(data['RETURN_NUMBER'] == 1)[0]
    goodReturns = numpy.zeros(data['X'].size, dtype=numpy.uint8)
    for i in firstReturnIndex:
        numReturns = data['NUMBER_OF_RETURNS'][i]
        if numReturns == 1:
            goodReturns[i] = 1
        if numReturns > 1:
            t1 = data['TIMESTAMP'][i]
            t = data['TIMESTAMP'][i:i+numReturns]
            returnNums = data['RETURN_NUMBER'][i:i+numReturns]
            if ((t1 * numReturns == numpy.sum(t)) and
                (numpy.sum(returnNums) == sum(range(1, numReturns+1)))):
                goodReturns[i:i+numReturns] = 1
    return data[goodReturns == 1]


def filterByPulse(points, density):
    """
    For each 10 x 10 m grid cell, the filtered number of pulses are selected
    randomly. The input 'density' is the desired number of pulses per cell
    (100m^2), so 300 = 3 pulses per m^2.
    """
    pixelSize = 10
    returnNumber = points['RETURN_NUMBER']
    numberOfReturns = points['NUMBER_OF_RETURNS']
    (x, y, z) = (points['X'], points['Y'], points['Z'])
    minX, maxX, minY, maxY, minZ, maxZ = rw_als_data.get_mmXYZ(x,y,z)
    if minX == int(minX):
        minX = minX - pixelSize
    if minY == int(minY):
        minY = minY - pixelSize
    tileXsize = numpy.ceil(maxX) - int(minX)
    tileYsize = numpy.ceil(maxY) - int(minY)
    nRows = int(numpy.ceil(tileYsize / pixelSize))
    nCols = int(numpy.ceil(tileXsize / pixelSize))
    
    firstReturnsToKeep = numpy.zeros(x.size, dtype=numpy.uint8)
    (row, col) = gridding_interp_fns.xyToRowCol(x, y, int(minX), int(numpy.ceil(maxY)), pixelSize)
    for r in range(nRows):
        for c in range(nCols):
            ind = numpy.where((row == r) & (col == c) & (returnNumber == 1))[0]
            numPulses = x[ind].size
            if numPulses > density:
                selection = numpy.random.choice(ind, size=density, replace=False)
                firstReturnsToKeep[selection] = 1
            else:
                firstReturnsToKeep[ind] = 1
    
    # Set pointsToKeep to 1 for all subsequent returns by using NUMBER_OF_RETURNS
    pointsToKeep = numpy.copy(firstReturnsToKeep)
    uniqueNums = numpy.unique(numberOfReturns[(numberOfReturns > 1) & (firstReturnsToKeep == 1)])
    for i in uniqueNums:
        for j in range(1, i):
            ind = numpy.where((numberOfReturns == i) & (firstReturnsToKeep == 1))[0]
            indsToKeep = ind + j
            indsToKeep = indsToKeep[indsToKeep < pointsToKeep.size]
            pointsToKeep[indsToKeep] = 1
    
    return points[pointsToKeep == 1]


def rec2las(outlaz, data, header, epsg):
    """
    Creates a laz file from a recArray of data, using the supplied epsg and
    header information.
    """
    las = laspy.create(point_format=header.point_format, file_version=header.version)
    las.header.scales = header.scales
    las.header.offset = header.offset
    las.header.add_crs(crs.CRS.from_user_input(epsg))
    las.return_num = data['RETURN_NUMBER']
    las.num_returns = data['NUMBER_OF_RETURNS']
    las.gps_time = data['TIMESTAMP']
    las.intensity = data['INTENSITY']
    las.classification = data['CLASSIFICATION'] 
    las.x = data['X']
    las.y = data['Y']
    las.z = data['Z']
    las.write(outlaz)


dirNames = [r"baradine_201407_lid1",
            r"goondiwindi_201501_lid1",
            r"gwabegar_201401_lid1",
            r"katoomba_201804_lid1"]

epsgCodes = [28355, 28356, 28355, 28356]

# This lists the mean pulse densities for LID1, LID2/Apron by survey, which was
# caclulated using get_pulse_densities.py and stored in 
# C:\Users\Public\Documents\lid1_lid2_comparison\pulse_densities.csv
pulseDensityMaxMin = [[155, 31],
                      [168, 41],
                      [241, 33],
                      [227, 36]]

pulseDensityList = []
for p in pulseDensityMaxMin:
    b = int((p[0] - p[1]) / 5.0)
    pList = list(range(p[0], p[1], -b))
    pList[-1] = p[1]
    pulseDensityList.append(pList[1:])

for i, dirName in enumerate(dirNames):
    inDir = os.path.join(r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data", dirName + r"\\laz")
    inLazList = glob.glob(os.path.join(inDir, "*.laz"))
    epsg = epsgCodes[i]

    outBase = os.path.join(r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\reduced_density_data", dirName)
    if not os.path.isdir(outBase):
        os.mkdir(outBase)
    
    for pulseDensity in pulseDensityList[i]:
        
        outDir = os.path.join(outBase, 'laz_%03d'%pulseDensity)
        if not os.path.isdir(outDir):
            os.mkdir(outDir)

        for inLaz in glob.glob(os.path.join(inDir, '*.laz')):
            outLaz = os.path.join(outDir, os.path.basename(inLaz))
            if os.path.exists(outLaz) is True:
                print(outLaz)

            else:
                # Read in inLaz
                (data, header) = rw_als_data.laspy2rec(inLaz)

                # Filter to desired number of pulses per 100m^2
                data = filterByPulse(data, pulseDensity)

                # Write LAZ file with reduced point density
                rec2las(outLaz, data, header, epsg)
                print(outLaz)
            
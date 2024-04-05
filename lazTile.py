#!/usr/bin/env python
"""

conda create -n wbt -c conda-forge whitebox laspy lazrs-python

The merging with wbt or laspy does not work!

"""

import os
import sys
import argparse
import glob
import numpy
import laspy
import math
import shutil
from whitebox.whitebox_tools import WhiteboxTools


def append_to_las(in_las, out_las):
    """
    Appends in_laz points to out_laz
    """
    las = laspy.read(in_las)
    with laspy.open(out_las, mode='a') as outlas:
        outlas.append_points(las.points)


def main(inFile, inDir, outDir, tileSize):
    """
    Main function to tile the data
    """

    t = glob.glob(os.path.join(inDir, '*.laz'))
    outLas = os.path.join(inDir, 'TEST.laz')
    
    shutil.copyfile(t[0], outLas)
    for inLas in t[1:]:
        append_to_las(inLas, outLas)
    
    sys.exit()
    
    # Get all input files
    if inFile is None:
        tempList = glob.glob(os.path.join(inDir, "*"))
        lasList = []
        for t in tempList:
            if t[-3:] in ['las', 'laz', 'LAS', 'LAZ']:
                lasList.append(t)
    else:
        lasList = [inFile]
        inDir = os.path.dirname(inFile)
    
    # Tile files
    #for lasFile in lasList:
    #    wbt = WhiteboxTools()
    #    wbt.set_verbose_mode(False)
    #    wbt.lidar_tile(lasFile, width=tileSize, height=tileSize, 
    #                   origin_x=0.0, origin_y=0.0, min_points=0)
    
    # Get all tile names
    tileList = []
    for lasName in lasList:
        subdir = lasName[:-4]
        tileList += glob.glob(os.path.join(subdir, '*.laz'))
    
    # Get the extents of each tile
    xyList = []
    for lasFile in tileList:
        las = laspy.read(lasFile)
        x_min = int(numpy.min(las.x) / tileSize) * tileSize
        y_min = int(numpy.min(las.y) / tileSize) * tileSize
        xy = 'x%i_y%i'%(x_min, y_min)
        xyList.append(xy)
        x_max = math.ceil(numpy.max(las.x) / tileSize) * tileSize
        y_max = math.ceil(numpy.max(las.y) / tileSize) * tileSize
        x_size = x_max - x_min
        y_size = y_max - y_min
        print(os.path.basename(lasFile), xy, x_min, y_min, x_size, y_size)
    tileList = numpy.array(tileList)
    xyList = numpy.array(xyList)
    
    # Merge or copy tiles into outDir
    #for xy in numpy.unique(xyList):
    #    outLas = os.path.join(outDir, '%s.laz'%xy)
    #    t = tileList[xyList == xy]
    #    if len(t) == 1:
    #        shutil.copyfile(t[0], outLas)
    #    else:
    #        shutil.copyfile(t[0], outLas)
    #        for inLas in t[1:]:
    #            append_to_las(inLas, outLas)
    
    # Check the extents of the final tiles
    for lasFile in glob.glob(os.path.join(outDir, "*.laz")):
        las = laspy.read(lasFile)
        x_min = int(numpy.min(las.x) / tileSize) * tileSize
        y_min = int(numpy.min(las.y) / tileSize) * tileSize
        x_max = math.ceil(numpy.max(las.x) / tileSize) * tileSize
        y_max = math.ceil(numpy.max(las.y) / tileSize) * tileSize
        x_size = x_max - x_min
        y_size = y_max - y_min
        print(os.path.basename(lasFile), x_min, y_min, x_size, y_size)


def getCmdargs():
    """
    Get the command line arguments.
    """
    p = argparse.ArgumentParser(
            description=("Tiles LAS/LAZ files into smaller files."))
    p.add_argument("-f", "--inFile", dest="inFile", default=None,
                   help=("Input LAS/LAZ file"))
    p.add_argument("-i", "--inDir", dest="inDir", default=None,
                   help=("Input directory with LAS/LAZ files"))
    p.add_argument("-o", "--outDir", dest="outDir", default=None,
                   help=("Output directory for LAZ tiles"))
    p.add_argument("-t", "--tileSize", dest="tileSize", default=200, type=int,
                   help=("Tile size (width and hieght) in metres (default=200)"))
    cmdargs = p.parse_args()
    if (cmdargs.inFile is None and cmdargs.inDir is None):
        p.print_help()
        print("Must name input file or directory and output file.")
        sys.exit()
    if cmdargs.outDir is None:
        p.print_help()
        print("Must name input file or directory and output file.")
        sys.exit()
    return cmdargs


if __name__ == "__main__":
    cmdargs = getCmdargs()
    main(cmdargs.inFile, cmdargs.inDir, cmdargs.outDir, cmdargs.tileSize)
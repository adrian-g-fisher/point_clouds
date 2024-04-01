#!/usr/bin/env python
"""

conda create -n wbt -c conda-forge whitebox

"""

import os
import sys
import argparse
import glob
from whitebox.whitebox_tools import WhiteboxTools


def main(inFile, inDir, outDir, tileSize):
    """
    Main function to tile the data
    """
    if inFile is None:
        tempList = glob.glob(os.path.join(inDir, "*"))
        lasList = []
        for t in tempList:
            if t[-3:] in ['las', 'laz', 'LAS', 'LAZ']:
                lasList.append(t)
    
    singleLaz = os.path.join(outDir, "merged.laz")
    wbt = WhiteboxTools()
    wbt.lidar_join(';'.join(lasList), singleLaz)

    # Tile files
    wbt = WhiteboxTools()
    wbt.lidar_tile(singleLaz, width=tileSize, height=tileSize, 
                   origin_x=0.0, origin_y=0.0, min_points=0)


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
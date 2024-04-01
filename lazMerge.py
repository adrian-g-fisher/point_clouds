#!/usr/bin/env python
"""
Merges multiple LAZ files and creates a new LAZ file with all data.

conda create -n als -c conda-forge laspy lazrs-python

"""

import os
import sys
import numpy
import argparse
import laspy
import glob


def laz2rec(inLaz):
    with laspy.open(inLaz) as fh:
        las = fh.read()
    header = las.header
    data = numpy.rec.fromarrays([las.return_num, las.num_returns, las.gps_time,
                                 las.intensity, las.classification, las.x,
                                 las.y, las.z],
                                 names=["RETURN_NUMBER", "NUMBER_OF_RETURNS",
                                        "TIMESTAMP", "INTENSITY",
                                        "CLASSIFICATION", "X", "Y", "Z"],
                                 formats=["u1", "u1", "<f8", "<i4", "u1", "<f8",
                                          "<f8", "<f8"])
    return data, header


def merge_structured_arrays(array1, array2):
    n1 = len(array1)
    n2 = len(array2)
    array_out = array1.copy()
    array_out.resize(n1 + n2)
    array_out[n1:] = array2
    return array_out


def rec2laz(outlaz, data, header):
    """
    Creates a laz file from a recArray of data, using the supplied epsg and
    header information.
    """
    las = laspy.create(point_format=header.point_format, file_version=header.version)
    las.header.scales = header.scales
    las.header.offset = header.offset
    las.return_num = data['RETURN_NUMBER']
    las.num_returns = data['NUMBER_OF_RETURNS']
    las.gps_time = data['TIMESTAMP']
    las.intensity = data['INTENSITY']
    las.classification = data['CLASSIFICATION'] 
    las.x = data['X']
    las.y = data['Y']
    las.z = data['Z']
    las.write(outlaz)


def merge_laz(inDir, outLaz):
    """
    Main function to merge the data
    """
    
    if os.path.exists(outLaz):
        print("ERROR: outLaz cannot exist")
        sys.exit()
    
    inLazList = glob.glob(os.path.join(inDir, "*.las"))
    for i, lazFile in enumerate(inLazList):
        if i == 0:
            outData, outHeader = laz2rec(lazFile)
        else:
            data, header = laz2rec(lazFile)
            outData = merge_structured_arrays(outData, data)
    
    rec2laz(outLaz, outData, outHeader)


def getCmdargs():
    """
    Get the command line arguments.
    """
    p = argparse.ArgumentParser(
            description=("Merges multiple LAZ files and creates a new LAZ file with all data."))
    p.add_argument("-i", "--inDir", dest="inDir", default=None,
                   help=("Input directory with LAZ files to merge"))
    p.add_argument("-o", "--outLaz", dest="outLaz", default=None,
                   help=("Output LAZ file"))
    cmdargs = p.parse_args()
    if (cmdargs.inDir is None or cmdargs.outLaz is None):
        p.print_help()
        print("Must name input directory and output LAZ file.")
        sys.exit()
    return cmdargs


if __name__ == "__main__":
    cmdargs = getCmdargs()
    merge_laz(cmdargs.inDir, cmdargs.outLaz)
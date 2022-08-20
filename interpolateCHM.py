#!/usr/bin/env python

"""
This interpolates a canopy height model (chm) from a maximum height (maxh)
image. It also uses a ground return image so that it can determine the extent of
the data. The method is as follows:
- nodata pixels are filled by using natural neigbour interpolation.
- canopy pits are identified using the adadptive mean filter (AMF) method
  modified from Shamsoddini et al. (2013).
- identified canopy pits are filled using natural neigbour interpolation. 
- any pixel further than 10 pixels from a ground return is considered beyond
  the data extent and is given the nodata value.

References

Shamsoddini, A., Turner, R., & Trinder, J.C. (2013). Improving lidar-based
forest structure mapping with crown-level pit removal. Journal of Spatial
Science, 58, 29-51, http://dx.doi.org/10.1080/14498596.2012.759092
"""

from __future__ import print_function, division

import sys
import os
import subprocess
import argparse
import numpy as np
from osgeo import gdal
from rios import applier
from scipy import ndimage
from pynninterp import NaturalNeighbour


def circular_element(s):
    """
    Makes a circular structuring element of width s.
    """
    se = np.ones((s, s,), dtype=np.uint8)
    se[int(s/2.0), int(s/2.0)] = 0
    se = ndimage.distance_transform_edt(se)
    se = (se <= int(s/2.0)).astype(np.uint8)
    return se


def make_chm(info, inputs, outputs, otherargs):
    """
    Called through RIOS to create the canopy height model. This assumes the
    nodata value of the maxh image is -999 and the nodata value of the ground
    point image is 255, which are the standard product values.
    
    It also ignores any max height pixels > 100 m above ground level (we assume
    they are bird hits).
    """
    # Interpolate the CHM if there are missing values
    maxh = inputs.maxh[0]
    maxh[maxh > 100] = -999
    if np.sum(maxh == -999) > 0:
        (xArray, yArray) = info.getBlockCoordArrays()
        xVals = xArray[maxh > -999]
        yVals = yArray[maxh > -999]
        zVals = maxh[maxh > -999].astype(np.float64)
        if zVals.size >= 3:
            chm = NaturalNeighbour(xVals, yVals, zVals, xArray, yArray)
            chm[np.isnan(chm)] = -999
            chm[chm < 0] = 0
            chm = chm.astype(np.float32)
        else:
            chm = maxh
    else:
        chm = maxh
    
    # Determine the extent of the data using the ground returns
    grd = inputs.grd[0]
    grd[grd == 255] = 0
    grd = ndimage.binary_closing(grd > 0, structure=circular_element(21))
    grd = ndimage.binary_dilation(grd > 0,
                                  structure=ndimage.generate_binary_structure(2, 2))
    grd = grd.astype(np.float32)
    grd[grd == 0] = -999
    chm[grd == -999] = -999
    
    # Output the chm with pits
    outputs.chmwithpits = np.array([chm]).astype(np.float32)
    
    # Identify and fill canopy pits
    if np.sum(chm != -999) > 0:
        
        # Calculate mean of of 8 neighbours while ignoring nodata values
        kernel = np.ones((3, 3), dtype=np.float32)
        kernel[1, 1] = 0
        chmValid = (chm > -999).astype(np.float32)
        chmValidCount = ndimage.convolve(chmValid, kernel)
        chmValidCount[chmValidCount == 0] = -1
        chmZerod = np.where(chm == -999, 0, chm)
        chmValidSum = ndimage.convolve(chmZerod, kernel)
        meanCHM = np.divide(chmValidSum, chmValidCount)
        meanCHM[chmValidCount == -1] = -999
        
        # Calculate stdev of 8 neighbours while ignoring nodata values
        chmValidSquaredSum = ndimage.convolve(np.power(chmZerod, 2), kernel)
        meanSquaredCHM = np.divide(chmValidSquaredSum, chmValidCount)
        diff = meanSquaredCHM - np.power(meanCHM, 2)
        diff[diff < 0] = 0
        stdevCHM = np.power(diff, 0.5)
        stdevCHM[stdevCHM == 0] = 0.001
        
        # Calculate the normalised value (nv) of the CHM
        # (the number of std devs that the CHM is different from its 8 neighbours)
        nv = (chm - meanCHM) / stdevCHM
        
        # New idea:
        # - simplify si as the mean of the nv 8 neighbours
        # - calculate gn as the count of 8 neighbours < 0.5 
        # - define pits as:
        #       lower nv than si
        #       negative nv
        #       count of 8 ground neighbours < 4
        #       mean of 8 neighbours > 0.5 m
        # - define bird strikes as:
        #       greater nv than si
        #       large positive nv
        #       chm > 20
        si = ndimage.convolve(nv, kernel) / np.sum(kernel)
        gn = ndimage.convolve((chm < 0.5).astype(np.uint8), kernel)
        pits = ((nv < 0) & (nv < si) & (gn < 5) & (meanCHM > 0.5)).astype(np.uint8)
        pits[(nv > 2) & (nv > si) & (chm > 50)] = 2 # bird strikes
        
        # Interpolate the CHM over canopy pits
        if np.sum(pits > 0) > 0:
            truePixels = np.where(maxh > -999, 1, 0)
            truePixels[grd == -999] = 0
            truePixels[pits > 0] = 0
            (xArray, yArray) = info.getBlockCoordArrays()
            xVals = xArray[truePixels == 1].astype(np.float64)
            yVals = yArray[truePixels == 1].astype(np.float64)
            zVals = maxh[truePixels == 1].astype(np.float64)
            if zVals.size >= 3:
                chm = NaturalNeighbour(xVals, yVals, zVals, xArray, yArray)
                chm[np.isnan(chm)] = -999
                chm[chm < 0] = 0
                chm = chm.astype(np.float32)
            else:
                chm = maxh
        
        # Set nodata pixels
        pits[grd == -999] = 0
        nv[grd == -999] = -999
        si[grd == -999] = -999
        chm[grd == -999] = -999
        meanCHM[grd == -999] = -999
        gn[grd == -999] = 0
    
    else:
        nv = np.zeros_like(chm) - 999
        si = np.zeros_like(chm) - 999
        meanCHM = np.zeros_like(chm) -999
        pits = np.zeros_like(chm)
        gn = np.zeros_like(chm)
    
    grd[grd == -999] = 0
    
    outputs.nv = np.array([nv]).astype(np.float32)
    outputs.si = np.array([si]).astype(np.float32)
    outputs.pits = np.array([pits]).astype(np.uint8)
    outputs.chm = np.array([chm]).astype(np.float32)
    outputs.mask = np.array([grd]).astype(np.uint8)
    outputs.mean = np.array([meanCHM]).astype(np.float32)
    outputs.gn = np.array([gn]).astype(np.uint8)


def run_maxh2chm(maxh_image, grd_image, chm_image):
    """
    This sets up RIOS to make the CHM.    
    """
    infiles = applier.FilenameAssociations()
    infiles.maxh = maxh_image
    infiles.grd = grd_image
    outfiles = applier.FilenameAssociations()
    outfiles.chm = chm_image
    outfiles.chmwithpits = chm_image.replace('chm', 'chmwithpits')
    outfiles.nv = chm_image.replace('chm', 'nv')
    outfiles.si = chm_image.replace('chm', 'si')
    outfiles.pits = chm_image.replace('chm', 'pits')
    outfiles.mask = chm_image.replace('chm', 'mask')
    outfiles.mean = chm_image.replace('chm', 'mean')
    outfiles.gn = chm_image.replace('chm', 'gn')
    controls = applier.ApplierControls()
    controls.setStatsIgnore(-999)
    controls.setOverlap(11)
    
    if chm_image[-3:] == "tif":
        controls.setOutputDriverName('GTiff')
    
    otherargs = applier.OtherInputs()
    applier.apply(make_chm, infiles, outfiles, otherargs, controls)


def getCmdargs():
    """
    Get the command line arguments.
    """
    p = argparse.ArgumentParser(
            description=("Creates a canopy height model from maximum height " +
                         "and ground return images"))
    p.add_argument("-m", "--maxh_image", dest="maxh_image", default=None,
                   help=("Maximum height image"))
    p.add_argument("-g", "--grd_image", dest="grd_image", default=None,
                   help=("Ground return image"))
    p.add_argument("-c", "--chm_image", dest="chm_image", default=None,
                   help=("Canopy height model image"))
    cmdargs = p.parse_args()
    if (cmdargs.maxh_image is None or
        cmdargs.grd_image is None or
        cmdargs.chm_image is None):
        p.print_help()
        print("Must name inputs and output")
        sys.exit()
    return cmdargs


if __name__ == "__main__":
    cmdargs = getCmdargs()
    run_maxh2chm(cmdargs.maxh_image, cmdargs.grd_image, cmdargs.chm_image)


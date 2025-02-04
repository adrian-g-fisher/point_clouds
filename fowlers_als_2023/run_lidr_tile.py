#!/usr/bin/env python

import os
import sys
import glob

for baseDir in glob.glob(os.path.join(r'D:\fowlers_als', '*')):
    inDir = os.path.join(baseDir, 'reprojected')
    outDir = os.path.join(baseDir, 'tiles')
    
    if os.path.exists(outDir) is False:
        os.mkdir(outDir)
    
    prefix = os.path.basename(glob.glob(os.path.join(inDir, '*'))[0]).replace('_z54.laz', '').replace('_', '')
    
    cmd = 'Rscript ../lidr_tile.r -i %s -o %s -t 1000 -x %s'%(inDir, outDir, prefix)
    os.system(cmd)
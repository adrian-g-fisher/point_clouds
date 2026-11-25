#!/usr/bin/env python

import os
import sys
import glob

for baseDir in glob.glob(os.path.join(r'D:\fowlers_lidar_comparison\drone', '*')):
    
    inDir = os.path.join(baseDir, 'terra_las')
    outDir = os.path.join(baseDir, 'tiles')
    
    if os.path.exists(outDir) is False:
        os.mkdir(outDir)
    
    prefix = os.path.basename(glob.glob(os.path.join(inDir, '*'))[0]).replace('.las', '')
    
    cmd = 'Rscript ../lidr_tile.r -i %s -o %s -t 250 -x %s'%(inDir, outDir, prefix)
    
    os.system(cmd)
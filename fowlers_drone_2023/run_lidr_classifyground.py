#!/usr/bin/env python

import os
import sys
import glob

for baseDir in glob.glob(os.path.join(r'D:\fowlers_lidar_comparison\drone', '*')):
    
    inDir = os.path.join(baseDir, 'tiles')
    outDir = os.path.join(baseDir, 'classified')
    
    if os.path.exists(outDir) is False:
        os.mkdir(outDir)
    
        origDir = os.path.join(baseDir, 'terra_las')
        prefix = os.path.basename(glob.glob(os.path.join(origDir, '*'))[0]).replace('.las', '')
    
        cmd = 'Rscript ../lidr_classifyground.r -i %s -o %s -x %s'%(inDir, outDir, prefix)
    
        os.system(cmd)
#!/usr/bin/env python
# First merge the B1 tiles together using lasmerge64

import os
import sys
import glob

baseDir = 'D:\\fowlers_als'
for inDir in glob.glob(os.path.join(baseDir, '*')):
    origDir = os.path.join(inDir, 'original')
    reprojDir = os.path.join(inDir, 'reprojected')
    if os.path.exists(reprojDir) is False:
        os.mkdir(reprojDir)
    for inLaz in glob.glob(os.path.join(origDir, '*.las')):
        outLaz = os.path.join(reprojDir, os.path.basename(inLaz).replace('.las', '_z54.laz'))
        if os.path.exists(outLaz) is False:
            cmd = r'C:\Users\z9803884\Documents\LAStools\bin\las2las -i %s -epsg 7855 -target_epsg 7854 -o %s'%(inLaz, outLaz)
            os.system(cmd)
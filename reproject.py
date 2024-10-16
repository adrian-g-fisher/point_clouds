#!/usr/bin/env python

import os
import sys
import glob


baseDir = 'F:\\fowlersgap_als_2023'

for inDir in glob.glob(os.path.join(baseDir, '*')):
    badDir = os.path.join(inDir, 'no_duplicates')
    newDir = os.path.join(inDir, 'tiles')
    if os.path.exists(newDir) is False:
        os.mkdir(newDir)
    for inLaz in glob.glob(os.path.join(badDir, '*.laz')):
        outLaz = os.path.join(newDir, os.path.basename(inLaz))
        cmd = 'Rscript lidr_reproject.r -i %s -o %s'%(inLaz, outLaz)
        os.system(cmd)

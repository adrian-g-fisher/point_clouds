#!/usr/bin/env python

import os
import sys
import glob


for baseDir in glob.glob('D:\\fowlers_als\\fg_als*'):
    inDir = os.path.join(baseDir, 'tiles')
    outDir = os.path.join(baseDir, 'no_duplicates')
    
    if os.path.exists(outDir) is False:
        os.mkdir(outDir)
    
    for inFile in glob.glob(os.path.join(inDir, '*.laz')):
        outFile = os.path.join(outDir, os.path.basename(inFile))
        cmd = 'Rscript lidr_removeduplicates.r -i %s -o %s'%(inFile, outFile)
        os.system(cmd)
        
sys.exit()

for baseDir in glob.glob('D:\\fowlers_als\\fg_als*'):
    inDir = os.path.join(baseDir, 'las')
    outDir = os.path.join(baseDir, 'tiles')
    cmd = 'Rscript lidr_tile.r -i %s -o %s -t 1000'%(inDir, outDir)
    os.system(cmd)
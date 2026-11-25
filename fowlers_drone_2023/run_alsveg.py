#!/usr/bin/env python

import os
import sys
import glob

for baseDir in glob.glob(os.path.join(r'D:\fowlers_lidar_comparison\drone', '*')):
    
    inDir = os.path.join(baseDir, 'classified')
    outDir = os.path.join(baseDir, 'indexed')
    
    if os.path.exists(outDir) is False:
        os.mkdir(outDir)
    
    laz_flist = os.path.join(inDir, 'laslist.txt')
    with open(laz_flist, 'w') as f:
        for lasfile in glob.glob(os.path.join(inDir, '*.laz')):
            f.write(os.path.basename(lasfile) + '\n')
    
    site = os.path.basename(baseDir)
    if site == 'D2_030_045':
        tile_s = 25
        binSize = 5
    else:
        tile_s = 100
        binSize = 25
    
    site = site.replace('_0', '')
    if site == 'D2_10090':
        site = 'D20090'
    
    cmd = (r'python C:\Users\z9803884\Documents\GitLab\als\src\alsveg\scripts\run_preprocess_lasf.py ' +
            '--indir %s --outdr %s --laz_flist laslist.txt --tile_s %s --out_tile_s %s '%(inDir, outDir, tile_s, tile_s) +
            '--epsg 28354 --no-useSpatialD --ss ap --ii r1 --pp dr --proj %s --year 2023 --binSize %s'%(site, binSize))
    #os.system(cmd)
    
    laz_flist = os.path.join(outDir, 'laslist.txt')
    with open(laz_flist, 'w') as f:
        for lasfile in glob.glob(os.path.join(outDir, '*.laz')):
            f.write(os.path.basename(lasfile) + '\n')
    
    cmd = (r'python C:\Users\z9803884\Documents\GitLab\als\src\alsveg\scripts\run_alsveg_workflow.py ' +
            '--indir %s --laz_flist laslist.txt --tile_s %s --epsg 28354 --psize 0.05'%(outDir, tile_s))
    #os.system(cmd)
    
    mosaicDir = os.path.join(baseDir, 'mosaics')
    
    if os.path.exists(mosaicDir) is False:
        os.mkdir(mosaicDir)
    
    cmd = (r'python C:\Users\z9803884\Documents\GitLab\als\src\alsveg\scripts\run_mergeTiles.py ' +
            '--indir %s --outdr %s --laz_flist laslist.txt --tile_s %s --psize 0.05 --chm_psize 0.2'%(outDir, mosaicDir, tile_s))
    os.system(cmd)
    
#!/usr/bin/env python

import os
import glob
import shutil

dirNames = [r"baradine_201406_apron",   r"baradine_201407_lid1",
            r"goondiwindi_201501_lid1", r"goondiwindi_201506_lid2",
            r"gwabegar_201401_lid1",    r"gwabegar_201406_apron",
            r"katoomba_201804_lid1",    r"katoomba_201804_lid2"]

epsg_codes = [28355, 28355, 28356, 28356, 28355, 28355, 28356, 28356]

for i, dirName in enumerate(dirNames):
    inDir = os.path.join(r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data", dirName + r"\\laz")
    epsg = epsg_codes[i]
    cmd = "run_alsworkflow --indir %s --tile_s 2000 --epsg %i --psize 1.0 --split_fpc"%(inDir, epsg)
    os.system(cmd)

    outDir = os.path.join(r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data", dirName)
    cmd = "run_mergeTiles --indir %s --outdr %s --lazlist laztilelist.txt --tile_s 2000 --psize 1.0"%(inDir, outDir)
    os.system(cmd)

    # Remove tile products
    for lazFile in glob.glob(os.path.join(inDir, '*.laz')):
        tempDir = lazFile.replace('.laz', '')
        if os.path.exists(tempDir):
            shutil.rmtree(tempDir)

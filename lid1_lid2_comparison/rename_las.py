#!/usr/bin/env python

import os
import sys

# Iterate through each project and run preprocess_als to create laz files with
# names that can be used by the processing code.

indirs = [r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\baradine_201407_lid1\\las",
          r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\baradine_201406_apron\\las",
          r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\goondiwindi_201501_lid1\\las",
          r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\goondiwindi_201506_lid2\\las",
          r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\gwabegar_201401_lid1\\las",
          r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\gwabegar_201406_apron\\las",
          r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\katoomba_201804_lid1\\las",
          r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data\\katoomba_201804_lid2\\las"]

instruments = ["l1", "l1", "l1", "l1", "l1", "l1", "l4", "l4"]
projects = ["BARAL1", "BARAAP", "GOONL1", "GOONL2", "GWABL1", "GWABAP", "KATAL1", "KATAL2"]
dates = [20140703, 20140609, 20150107, 20150602, 20140116, 20140609, 20180412, 20180411]
epsg_codes = [28355, 28355, 28356, 28356, 28355, 28355, 28356, 28356]

for i, indir in enumerate(indirs):
    outdir = os.path.join(os.path.dirname(indir), "laz")
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    ii = instruments[i]
    proj = projects[i]
    date = dates[i]
    epsg = epsg_codes[i]
    cmd = ("preprocess_als --indir %s --outdr %s "%(indir, outdir) +
           "--tile_s 2000 --ntile_s 2000 --epsg %s --noCreateBA3 "%epsg +
           "--ss ap --ii %s --pp dr --proj %s --year %i"%(ii, proj, date))
    os.system(cmd)
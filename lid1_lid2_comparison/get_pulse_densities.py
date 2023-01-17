#!/usr/bin/env python

import os
import glob
import numpy
from rios import applier


def calcMeanPulseDensity(info, inputs, outputs, otherargs):
    pulseDensityPixels = inputs.pulseDensity[0][inputs.surveyAreas[0] == 1]
    
    if otherargs.values is None:
        otherargs.values = pulseDensityPixels
    else:
        otherargs.values = numpy.concatenate([otherargs.values, pulseDensityPixels])


surveyPolygons = r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\lidar_overlap.shp"

dirNames = [r"baradine_201406_apron",   r"baradine_201407_lid1",
            r"goondiwindi_201501_lid1", r"goondiwindi_201506_lid2",
            r"gwabegar_201401_lid1",    r"gwabegar_201406_apron",
            r"katoomba_201804_lid1",    r"katoomba_201804_lid2"]

with open('pulse_densities.csv', 'w') as f:
    f.write('image,pulse_density\n')

for i, dirName in enumerate(dirNames):
    inDir = os.path.join(r"C:\\Users\\Public\\Documents\\lid1_lid2_comparison\\input_data", dirName)
    pulseDensityImage = glob.glob(os.path.join(inDir, "*bb5*.tif"))[0]
    infiles = applier.FilenameAssociations()
    infiles.pulseDensity = pulseDensityImage
    infiles.surveyAreas = surveyPolygons
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    controls.setReferenceImage(pulseDensityImage)
    controls.setFootprintType(applier.BOUNDS_FROM_REFERENCE)
    otherargs.values = None
    applier.apply(calcMeanPulseDensity, infiles, outfiles, otherArgs=otherargs, controls=controls)
    with open('pulse_densities.csv', 'a') as f:
        f.write('%s,%f\n'%(os.path.basename(pulseDensityImage), numpy.mean(otherargs.values)))
# Point Cloud Processing
This is a collection of scripts for processing point clouds (LAS format files) from airborne lidar or drone photogrammetric surveys. Most of the are command line scripts that allow inputs and outputs to be specified. The scripts depend on several python packages, including RIOS and pylidar. These can be installed into a new Anaconda environment, with the following command:

`conda config --add channels conda-forge`<br />
`conda config --add channels rios`<br />
`conda create -n lidar rios scipy numba rios::pynninterp`<br />
`conda activate lidar`<br />
`pip install laspy[laszip]`

## Drone canopy height models
To create a digital surface model (DSM), digital elevation model (DEM), and canopy height model (CHM) from a LAZ file point cloud produced by structure from
motion processing of overlapping drone imagery (e.g. output from Pix4DMapper):

`droneLasProcessing.py --inLaz test.laz --epsg 32755 --projectName test --outDir ./ --pixelsize 0.05 --windowsize 40`

## Lidar canopy height models
To create canopy height models (CHMs) from images of maximum height above ground derived from airborne lidar surveys:

`interpolateCHM.py --maxh_image test_maxh.tif --grd_image test_grd.tif --chm_image test_chm.tif`

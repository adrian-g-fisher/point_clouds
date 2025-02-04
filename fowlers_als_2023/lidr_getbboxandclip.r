# conda create -n lidr -c conda-forge r-lidr r-future r-rgdal r-codetools r-optparse r-tidyverse

library(lidR)
library(sf)
suppressPackageStartupMessages(library(tidyverse))

# Get bbox, check in QGIS and edit CSV file
base_dir <- "D:/fowlers_als/fgpb01/original/"
filelist <- list.files(base_dir)
for (infile in filelist) {
	print(infile)
	inlas <- paste0(base_dir, infile)
	las <- readLAS(inlas)
	print(st_bbox(las))
	
	#las <- readLAScatalog(inlas)
	#ctg <- catalog_boundaries(las, concavity = 1, length_threshold = 15)
	#ctg_sf <- st_as_sf(ctg)
	#outshape <- str_replace(inlas, ".las", "_boundary.shp")
	#st_write(ctg_sf, outshape)
}

# Clip to new bbox


#bbox_table <- read.csv('bbox.csv')
#base_dir <- 'D:/fowlers_als/fgpb01/reprojected/'
#for (infile in bbox_table$file) {
#	values = bbox_table[bbox_table$file == infile,]
#   inlas <- paste0(base_dir, infile)
#	las <- readLAS(inlas)
#	subset <- clip_rectangle(las, values$new_xmin, values$new_ymin,
#								  values$new_xmax, values$new_ymax)
#	outfile = str_replace(inlas, '.laz', '_clip.laz')
#	writeLAS(subset, file=outfile)
#}
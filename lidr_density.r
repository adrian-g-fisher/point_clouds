# conda create -n lidr -c conda-forge r-lidr r-future r-rgdal r-codetools r-optparse

library(lidR)
library(terra)
library(optparse)
library(future)

option_list = list(
	make_option(c("-i", "--inDir"), type="character", default=NULL, 
                help="Directory with input las files"),
	make_option(c("-o", "--outFile"), type="character", default=NULL, 
                help="Output density TIF file"),
	make_option(c("-p", "--pixelSize"), type="numeric", default=0.1, 
                help="Pixel size (m) for output TIF file"))	

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$inDir)){
  print_help(opt_parser)
  stop("Must supply inDir")}
if (is.null(opt$outFile)){
  print_help(opt_parser)
  stop("Must supply outFile")}

inDir <- opt$inDir
outFile <- opt$outFile
pixelSize <- opt$pixelSize

points <- readLAScatalog(inDir)
plan(multisession, workers=8L) # 8 cores
outgrid <- grid_density(points, res=pixelSize)
writeRaster(outgrid, outFile)
# conda create -n lidr -c conda-forge r-lidr

library(lidR)
library(optparse)

option_list = list(
	make_option(c("-i", "--inDir"), type="character", default=NULL, 
                help="Directory with input las files"),
	make_option(c("-o", "--outFile"), type="character", default=NULL, 
                help="Output mosaic TIF file"),
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

RGBZ <- function(r, g, b, z) {
	maxZ <- which.max(z)
	bands = list(R = r[maxZ], G = g[maxZ], B = b[maxZ])
	return(bands)}

points <- readLAScatalog(inDir)
RGB <- grid_metrics(points, ~RGBZ(R, G, B, Z), res=pixelSize)
outFile <- as.raster(RGB)
# conda create -n lidr -c conda-forge r-lidr

library(lidR)
library(optparse)

option_list = list(
	make_option(c("-i", "--inDir"), type="character", default=NULL, 
                help="Directory with input las files"),
	make_option(c("-o", "--outFile"), type="character", default=NULL, 
                help="Output mosaic TIF file"),
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

# returns the RGB values for the highest points
RGBZ <- function(r,g,b,z) {
  bands = list(
    R = r[which.max(z)],
    G = g[which.max(z)],
    B = b[which.max(z)]
  )
  return(bands)
}

points <- readLAScatalog(inDir)
ortho <- grid_metrics(points, ~RGBZ(R,G,B,Z), res = 0.1)
writeRaster(ortho, outFile)

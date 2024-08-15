# conda create -n lidr -c conda-forge r-lidr r-future r-rgdal r-codetools r-optparse r-dplyr

library(lidR)
library(sf)
library(optparse)


option_list = list(
	make_option(c("-i", "--inFile"), type="character", default=NULL, 
                help="Input las/laz file"),
	make_option(c("-o", "--outFile"), type="character", default=NULL, 
                help="Output las/laz file"),
    make_option(c("-e", "--inEPSG"), type="numeric", default=7855, 
                help="Input file EPSG code"),
    make_option(c("-f", "--outEPSG"), type="numeric", default=7854, 
                help="Output file EPSG code"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$inFile)){
  print_help(opt_parser)
  stop("Must supply inFile")}
if (is.null(opt$outFile)){
  print_help(opt_parser)
  stop("Must supply outFile")}

inFile <- opt$inFile
outFile <- opt$outFile
inEPSG <- opt$inEPSG
outEPSG <- opt$outEPSG

las <- readLAS(inFile)
epsg(las) <- inEPSG
las2 <- st_transform(las, st_crs(outEPSG))
writeLAS(las2, outFile)
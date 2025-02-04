# conda create -n lidr -c conda-forge r-lidr r-future r-rgdal r-codetools r-optparse
# 
# If the input files have an overlap they might need to be clipped to non-overlapping areas first
# > lasclipRectangle(ctg, xleft, ybottom, xright, ytop, ofile = "")

library(lidR)
library(optparse)
library(future)

option_list = list(
	make_option(c("-i", "--inDir"), type="character", default=NULL, 
                help="Directory with input las files"),
	make_option(c("-o", "--outDir"), type="character", default=NULL, 
                help="Directory for output las files"),
	make_option(c("-t", "--tileSize"), type="numeric", default=NULL, 
                help="Tile size (m) for output las files"),
	make_option(c("-x", "--prefix"), type="character", default=NULL, 
                help="Prefix for output las files"))
				
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$inDir)){
  print_help(opt_parser)
  stop("Must supply inDir")}
if (is.null(opt$outDir)){
  print_help(opt_parser)
  stop("Must supply outDir")}
if (is.null(opt$tileSize)){
  print_help(opt_parser)
  stop("Must supply tileSize")}
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Must supply prefix")}
  
inDir <- opt$inDir
outDir <- opt$outDir
tileSize <- opt$tileSize
prefix <- opt$prefix

ctg <- readLAScatalog(inDir)
plan(multisession, workers=8L) # 8 cores
opt_chunk_buffer(ctg) <- 0
opt_chunk_size(ctg) <- tileSize
opt_laz_compression(ctg) <- TRUE
opt_chunk_alignment(ctg) <- c(0, 0)
opt_output_files(ctg) <- paste0(outDir, "/", prefix, "_{XLEFT}_{YTOP}")
newctg <- catalog_retile(ctg)
# conda create -n lidr -c conda-forge r-lidr r-future r-rgdal r-codetools r-optparse

library(lidR)
library(optparse)
library(future)

option_list = list(
	make_option(c("-i", "--inDir"), type="character", default=NULL, 
                help="Directory with input las files"),
	make_option(c("-o", "--outDir"), type="character", default=NULL, 
                help="Directory for output las files"))
				
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$inDir)){
  print_help(opt_parser)
  stop("Must supply inDir")}
if (is.null(opt$outDir)){
  print_help(opt_parser)
  stop("Must supply outDir")}

inDir <- opt$inDir
outDir <- opt$outDir

ctg = readLAScatalog(inDir)
plan(multisession, workers=8L) # 8 cores
opt_chunk_buffer(ctg) <- 20
opt_laz_compression(ctg) <- TRUE
opt_output_files(ctg) <- paste0(outDir, "/{*}_classified")
ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out=length(ws))
ctg <- classify_ground(ctg, algorithm=pmf(ws=ws, th=th))
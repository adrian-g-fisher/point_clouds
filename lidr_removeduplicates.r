library(lidR)
library(optparse)
library(future)

option_list = list(
	make_option(c("-i", "--inFile"), type="character", default=NULL, 
                help="Input laz file with duplicates"),
	make_option(c("-o", "--outFile"), type="character", default=NULL, 
                help="Output laz file without duplicates"))
				
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

las <- readLAS(inFile)
newlas <- filter_duplicates(las)
writeLAS(newlas, outFile)
# conda create -n lidr -c conda-forge r-lidr
# 
# If the input files have an overlap (duplicate points), they might be flagged,
# in which case they can be filtered out with:
# > opt_filter(ctg) <- "-drop_withheld"
#
# Otherwise, the input files might need to be clipped to non-overlapping areas.
# > lasclipRectangle(ctg, xleft, ybottom, xright, ytop, ofile = "")
#
# Or we can remove duplicates using a filter
# > filter_duplicates(ctg)

library(lidR)
library(optparse)

option_list = list(
	make_option(c("-i", "--inDir"), type="character", default=NULL, 
                help="Directory with input las files"),
	make_option(c("-o", "--outDir"), type="character", default=NULL, 
                help="Directory for output las files"),
	make_option(c("-t", "--tileSize"), type="numeric", default=NULL, 
                help="Tile size (m) for output las files"),
	make_option(c("-c", "--classify"), action="store_true", default=FALSE, 
                help="Classify ground with progressive morphological filter"))
				
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

inDir <- opt$inDir
outDir <- opt$outDir
tileSize <- opt$tileSize
classify <- opt$classify

ctg = readLAScatalog(inDir)
opt_chunk_size(ctg) <- tileSize
opt_laz_compression(ctg) <- TRUE
opt_filter(ctg) <- ""
opt_chunk_alignment(ctg) <- c(0, 0)
opt_output_files(ctg) <- paste0(outDir, "/retile_{XLEFT}_{YBOTTOM}")

if (classify == TRUE) {
	opt_chunk_buffer(ctg) <- 20
	ws <- seq(3, 12, 3)
	th <- seq(0.1, 1.5, length.out=length(ws))
	ctg <- classify_ground(ctg, algorithm=pmf(ws=ws, th=th))}

opt_chunk_buffer(ctg) <- 0
newctg = catalog_retile(ctg)

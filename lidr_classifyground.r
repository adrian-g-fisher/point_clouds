# conda create -n lidr -c conda-forge r-lidr r-future r-rgdal r-codetools r-optparse r-dplyr

library(lidR)
library(optparse)
library(future)
suppressPackageStartupMessages(library(dplyr))


lowest <- function(res=1) {
  stopifnot(is.numeric(res), length(res) == 1L, res > 0)
  
  f <- function(las) {
    r <- grid_metrics(las, ~.I[which.min(Z)], res=res)
    return(na.omit(r[]))
  }
  
  class(f) <- lidR:::LIDRALGORITHMDEC # Using ::: because not public yet
  return(f)
}


minZ_pmf <- function(las, res=0.1) {
  if (is(las, "LAS")) {
    minZ <- decimate_points(las, lowest(res))
	ws <- seq(3, 12, 3)
	th <- seq(0.1, 1.5, length.out=length(ws))
	minZ <- classify_ground(minZ, algorithm=pmf(ws=ws, th=th))
	joined <- las@data %>% left_join(minZ@data, by=c('X','Y','ReturnNumber'))
	joined[is.na(joined)] <- 0
	las$Classification <- as.integer(joined$Classification.y)
	return(las)
  }
  
  if (is(las, "LAScatalog")) {
    options <- list(need_output_file = TRUE, need_buffer = TRUE)
    las <- catalog_map(las, minZ_pmf, res=res)
    return(las)
  }
}


option_list = list(
	make_option(c("-i", "--inDir"), type="character", default=NULL, 
                help="Directory with input las files"),
	make_option(c("-o", "--outDir"), type="character", default=NULL, 
                help="Directory for output las files"),
	make_option(c("-r", "--res"), type="numeric", default=0.1, 
                help="Resolution (m) for minimum elevation grid (default=0.1)"))

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
res <- opt$res

ctg = readLAScatalog(inDir)
plan(multisession, workers=8)
opt_chunk_buffer(ctg) <- 20
opt_laz_compression(ctg) <- TRUE
opt_output_files(ctg) <- paste0(outDir, "/{*}_classified")
output <- minZ_pmf(ctg, res=res)
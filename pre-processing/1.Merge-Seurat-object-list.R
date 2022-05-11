#!/usr/bin/env Rscript

# Merge list of individual Seurat objects, or read individual object

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Comma separated string (no spaces) with multiple paths to input Seurat objects [or] Path to list of objects'
  ), 
  make_option(
    c("-l", "--list"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Wether the input_seurat is a list of Seurat objects'
  ), 
  make_option(
    c("-o", "--output_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to output Seurat object'
  )
  )
  
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Seurat))
# if input is list
if(opt$list == T){
  # 1. Load list of seurat objects
  seu_list <- readRDS(opt$input_seurat)
  # 2.Merge
  aggr <- Reduce(merge, seu_list)
  # if input is string of paths
} else {
  # 0. Edit input seurat object paths 
  paths <-  gsub(opt$input_seurat, pattern = " ", replacement = "")
  # 0. Split by comma
  paths <- unlist(strsplit( paths , ","))
  # 1. Read objects
  seu_list <- lapply(paths, readRDS)
  # 2. Merge
  aggr <- Reduce(merge, seu_list)
}
# 3. Save
saveRDS(aggr, opt$output_seurat)
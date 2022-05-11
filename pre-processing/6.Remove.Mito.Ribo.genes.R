#!/usr/bin/env Rscript

## Remove mitochondrial and ribosomal features from Seurat object or object list.

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input Seurat object'
  ), 
  make_option(
    c("-l", "--list"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Logical indicating if input object is a list of Seurar objects.'
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

## 1. Remove mitochondrial and ribosomal genes from cells
remove.ribo.mito.genes <- function(seu){
  ## Filter Cells based on their mito and ribo content
  
  print(paste0("Dimensions PRE-filtering: ", dim(seu)))
  
  seu <- seu[!grepl("^MT-", rownames(seu)) & !grepl("^RPS", rownames(seu)) & !grepl("^RPL", rownames(seu)), ]
  
  print(paste0("Dimensions POST-filtering: ", dim(seu)))
  
  seu
}
# 1. Load object
seu <- readRDS(opt$input_seurat)
if(opt$list == T){
  # 1. Filter mito ribo genes
  seu <- lapply(seu, function(seu_obj) remove.ribo.mito.genes(seu_obj))
}else{
  # 1. Filter mito ribo genes
  seu <- remove.ribo.mito.genes(seu)
}
# 2. Save pre-processed object
saveRDS(seu, file = opt$output_seurat)

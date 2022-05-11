#!/usr/bin/env Rscript

# Quick pre-processing of a Seurat object. 

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
    help = 'Logical indicating if input object is a list of Seurat objects.'
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


seurat_preprocess <- function(seu){
  # Normalization
  seu <- SCTransform(seu, assay = "RNA", variable.features.n = 3000) # variable.features.n set to default
  
  # Scale Data (including all features)
  seu <- ScaleData(seu, assay = "SCT", features=NULL, do.scale = TRUE, do.center = TRUE)
  
  # Run PCA
  seu <- RunPCA(seu, assay = "SCT", verbose = F)
  # Neighbours 
  seu <- FindNeighbors(seu, reduction = "pca")
  # UMAP
  seu <- RunUMAP(seu, dims = 1:50)
  # Clusters 
  seu <- FindClusters(seu, resolution = c(seq(0.1, 1.1, 0.2), 1.3, 1.5, 2.0), alg = 1, graph.name = "SCT_snn")
  seu
}

# 1. Load object
seu <- readRDS(opt$input_seurat)

if(opt$list == T){
  # 2. Pre-process object
  seu <- lapply(seu, function(seu_obj)  seurat_preprocess(seu_obj))
}else{
  # 2. Pre-process object
  seu <- seurat_preprocess(seu)
}

# 3. Save pre-processed object
saveRDS(seu,opt$output_seurat)

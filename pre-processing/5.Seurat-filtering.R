#!/usr/bin/env Rscript

## Perform filtering of Seurat object or object list in terms of: nCount, nFeature, percent.mito, percent.ribo

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
    c("-m", "--mito_thres"),
    action = "store",
    default = 10,
    type = 'numeric',
    help = 'Maximum threshold of fraction of mitochondrial genes considered.'
  ),
  make_option(
    c("-r", "--ribo_thres"),
    action = "store",
    default = 10,
    type = 'numeric',
    help = 'Maximum threshold of fraction of ribosomal genes considered.'
  ),
  make_option(
    c("-a", "--nFeat_min"),
    action = "store",
    default = 800,
    type = 'numeric',
    help = 'Minimum number of Features for a cell to be considered.'
  ),
  make_option(
    c("-b", "--nCount_min"),
    action = "store",
    default = 1000,
    type = 'numeric',
    help = 'Min number of Counts for a cell to be considered.'
  ),
  make_option(
    c("-c", "--nCount_max"),
    action = "store",
    default = 50000,
    type = 'numeric',
    help = 'Max number of Counts for a cell to be considered.'
  ),
  make_option(
    c("-o", "--output_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to output Seurat object, or object list.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Seurat))


## 1. Filter cells based on mitochondrial and ribosomal gene content
subset.cells <- function(seu, thr.mt, thr.rp, min.thr.nFeat.RNA, min.thr.nCount.RNA, max.thr.nCount.RNA){
  ## Filter Cells based on their mito and ribo content
  
  # Checks
  if(!(all(c("percent.mt.RNA", "percent.rp.RNA", "nFeature_RNA", "nCount_RNA") %in% names(seu@meta.data)))) stop ("Not all parameters present in Seurat object meta.data")
  ## Replace NAs with a value, if not they get filtered out
  seu@meta.data$percent.mt.RNA [ is.na(seu@meta.data$percent.mt.RNA) ] <- 0
  seu@meta.data$percent.rp.RNA [ is.na(seu@meta.data$percent.rp.RNA) ] <- 0
  
  print(paste0("Dimensions PRE-filtering: ", dim(seu)))
  
  seu <- subset(seu, percent.mt.RNA < thr.mt & percent.rp.RNA < thr.rp &
                  nFeature_RNA > min.thr.nFeat.RNA & nCount_RNA > min.thr.nCount.RNA & nCount_RNA < max.thr.nCount.RNA)
  
  print(paste0("Dimensions POST-filtering: ", dim(seu)))
  
  seu
}

# 0. params
mito_thres = opt$mito_thres
ribo_thres = opt$ribo_thres
nFeat_min = opt$nFeat_min
nCount_min = opt$nCount_min
nCount_max = opt$nCount_max

# 1. Load object
seu <- readRDS(opt$input_seurat)
if(opt$list == T){
  for (obj in seu) {
    # calc rna, and mt %
    obj[["percent.mt.RNA"]] <- PercentageFeatureSet(obj, assay="RNA",  pattern = "^MT-")
    obj[["percent.rp.RNA"]] <- PercentageFeatureSet(obj, assay = "RNA", pattern = "^RPS") + PercentageFeatureSet(obj, assay = "RNA", pattern = "^RPL")
  }
  
  # 1. Filter cells
  seu <- lapply(seu, function(seu_obj) 
    subset.cells(seu_obj, 
                 thr.mt = mito_thres, 
                 thr.rp = ribo_thres, 
                 min.thr.nFeat.RNA = nFeat_min, 
                 min.thr.nCount.RNA = nCount_min, 
                 max.thr.nCount.RNA = nCount_max))
}else{
  # calc rna and mt %
  seu[["percent.mt.RNA"]] <- PercentageFeatureSet(seu, assay="RNA",  pattern = "^MT-")
  seu[["percent.rp.RNA"]] <- PercentageFeatureSet(seu, assay = "RNA", pattern = "^RPS") + PercentageFeatureSet(seu, assay = "RNA", pattern = "^RPL")
  
  
  # 1. Filter cells
  seu <-  subset.cells(seu, 
                       thr.mt = mito_thres, 
                       thr.rp = ribo_thres, 
                       min.thr.nFeat.RNA = nFeat_min, 
                       min.thr.nCount.RNA = nCount_min, 
                       max.thr.nCount.RNA = nCount_max)
}
# 3. Save pre-processed object
saveRDS(seu, file = opt$output_seurat)

#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(velocyto.R))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(scales))

option_list = list(
  make_option(
    c("-i", "--input_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input Seurat object, or list of Seurat objects'
  ),
  make_option(
    c("-v", "--input_velo"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input velocity estimate object, or list of objects.'
  ),
  make_option(
    c("-l", "--list"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Wether the input rds object is a list of Seurat objects or a single one.'
  ),
  make_option(
    c("-c", "--cl_annot"),
    action = "store",
    default = "SCT_snn_res.0.1",
    type = 'character',
    help = 'Cell cluster annotation to be considered in the UMAP color respresentation.'
  ),
  make_option(
    c("-o", "--output_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to save the RNA Velocity projection over embedding plots.'
  ), 
  make_option(
    c("-s", "--sample_name"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Seurat object sample name, provide if input is NOT a list.'
  ),
  make_option(
    c("-r", "--reduction"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Seurat object reduction to use for plotting.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(velocyto.R))

# Function #
color.vector <- function(seu, cl_annot){
  ## Obtain vector for cell coloring
  suppressPackageStartupMessages(require(scales, lib.loc = "/usr/local/lib/R/site-library/scales"))
  
  labels <- seu@meta.data[[cl_annot]]
  # print('cl_annot')
  # print(cl_annot)
  # print('labels')
  # print(labels)
  ncol <- length(unique(labels))
  x <- as.factor(labels)
  levels(x) <- 1:length(levels(x))
  x <- as.numeric(x)
  colvec <- hue_pal()(ncol)[x]
  
  # colvec <- hue_pal()(max(as.numeric(seu@meta.data[[cl_annot]]), na.rm = T))[as.numeric(seu@meta.data[[cl_annot]])]
  
  names(colvec) <- colnames(seu)
  return(colvec)
}

# 8. Velocity map on embedding
velocity.on.embedding.plot <- function(seu, reduction, pdf_path, grid_vals, velo.est, col_vec, sample_name, cl_annot){
  ## Generate plot of velocity on embedding
  #reduction=reduction # 'umap', 'umap_harmony_t.0.1'
  
  seu.umap = Embeddings(seu, reduction)
  pdf(pdf_path)
  for(val in grid_vals){
    show.velocity.on.embedding.cor(seu.umap, 
                                   velo.est,
                                   n=300,
                                   scale="sqrt",
                                   cell.colors = ac(col_vec, alpha=0.5),
                                   cex=0.8,
                                   arrow.scale=5,
                                   show.grid.flow=TRUE,
                                   min.grid.cell.mass=0.5,
                                   grid.n=val, # determined number of arrows!
                                   arrow.lwd=1,
                                   do.par=F,
                                   n.cores = 8,
                                   cell.border.alpha = 0.1,
                                   main=paste0("Sample: ",  sample_name, " Projection: UMAP ", cl_annot, " Grid n = ", val))
  }
  dev.off()
}



### Load files: 
# 1 Seurat object
print('loading')
print(opt$input_seurat)
seu = readRDS(opt$input_seurat)
print(names(seu@reductions))

# 2. Velocity estimates
print('loading')
print(opt$input_velo)
velo.est = readRDS(opt$input_velo)

# if list:
if(opt$list ==TRUE){
  #  iterate
  for (i in  1:length(seu)){
    seu_name = names(seu)[i]
    seu_obj = seu[[seu_name]]
    
    velo_obj = velo.est[[seu_name]]
    

    ## Subset Seurat object to velocity estimates object dimensions
    seu_obj <- seu_obj[velo_obj[["conv.nmat.norm"]]@Dimnames[[1]], velo_obj[["conv.nmat.norm"]]@Dimnames[[2]]]
    # 7. Color vector
    cl_annot = opt$cl_annot # colouring variable
    col.vec = color.vector(seu = seu_obj, cl_annot = cl_annot)
    # 8. Velocity map on embedding
    velocity.on.embedding.plot(seu = seu_obj, 
                                 pdf_path = paste0(opt$output_path, seu_name, "-Velocity-Plot-", cl_annot, ".pdf"),
                                 grid_vals = c(20, 40), 
                                 velo.est = velo_obj, col_vec = col.vec,
                                 sample_name = seu_name, cl_annot = cl_annot)
  }
  
  # if not list
  } else {
    
    print('not subsetting seu object')
    ## Subset Seurat object to velocity estimates object dimensions
    #seu <- seu[velo.est[["conv.nmat.norm"]]@Dimnames[[1]], velo.est[["conv.nmat.norm"]]@Dimnames[[2]]]
    
    # 7. Color vector
    print('getting col.vec')
    cl_annot = opt$cl_annot # colouring variable
    col.vec = color.vector(seu = seu, cl_annot = cl_annot)
    
    print('plotting...')
    # 8. Velocity map on embedding
    velocity.on.embedding.plot(seu = seu, 
                               reduction = opt$reduction,
                               pdf_path = paste0(opt$output_path, opt$sample_name, "-Velocity-Plot-", cl_annot, ".pdf"),
                               grid_vals = c(20, 40, 60), 
                               velo.est = velo.est, col_vec = col.vec,
                               sample_name = opt$sample_name, cl_annot = cl_annot)
    
  }
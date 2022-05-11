#!/usr/bin/env Rscript
rm(list=ls())

# Velocyto Output processing
## Generates RNA velocity representation in object embedding FOR CELL TYPE OBJECTS
## The difference here is that loom matrices have to be subseted to those cells composing the cell type object

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
    c("-t", "--cr_paths"),
    action = "store",
    default = "Velocyto_paths_CR_samples.csv",
    type = 'character',
    help = 'Path to the txt file containing CellRanger output directories. This file must contain two columns: `sample` indicating the object name and `cr_path` referring to the CellRanger output dir, where the velocyto dir is located.'
  ),
  make_option(
    c("-l", "--list"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Logical indicating if input object is a list of Seurat objects or an individual object.'
  ),
  make_option(
    c("-c", "--annots"),
    action = "store",
    default = "SCT_snn_res.0.1",
    type = 'character',
    help = 'Colour cells according to "annot" in UMAP RNA velocity plot. Character vector, string separated.'
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
    c('-z', "--subset_celltype"),
    action="store",
    default='All',
    type='character',
    help = 'the celltype to subset to prior to calculation. "All" if want all'
  ),
  make_option(
    c('-x', "--subset_disease_stage"),
    action='store',
    default='All',
    type='character',
    help='the disease stage to subset to prior to calculation. "All" if want all'
  ),
  make_option(
    c('-b', '--barcode_file_path'),
    action='store', 
    default=NA,
    type='character',
    help='path to file with barcodes which MUST be included in sampling'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# can't pass args with spaces
if (opt[['subset_disease_stage']]=='endStage') {
  opt[['subset_disease_stage']] <- 'end stage'
}
  

dir.create(opt$output_path, recursive = T)


### Functions ###

# 1. Read loom file
read.loom.file <- function(cr.dir){
  ## Read loom file from provided CellRanger output directory
  velo.dir = paste0(cr.dir, "velocyto/")
  loom.file = list.files(velo.dir, pattern =".loom", full.names=F)
  loom.matrices = read.loom.matrices(paste0(velo.dir, loom.file))
  loom.matrices
}

# 2. Rename loom files
rename.col.inner = function(x, lib.suffix){ 
  ## Rename loom matrix colnames to match Seurat object cell names
  new.names=rep(NA, length(x))
  for (i in 1:length(x)){
    # print('here')
    # print(x[i])
    sample=strsplit(x[i], '[:]')[[1]][1]
    barcode=strsplit(x[i], '[:]')[[1]][2]
    barcode=strsplit(barcode, '[x]')[[1]][1]
    # add -1 library preffix present in Seurat object cell names
    new.names[i] = toupper(paste0(barcode, lib.suffix))
  }
  new.names
}


rename.col.outer = function(y, lib.suffix){
  ## Rename loom matrices colnames to match Seurat object cell names
  new.col.names = rename.col.inner(colnames(y$spliced), lib.suffix)
  colnames(y$spliced) = new.col.names
  colnames(y$unspliced) = new.col.names
  colnames(y$ambiguous) = new.col.names
  y
}


# 3. Order and subset loom.matrices
order.loom.cells <- function(seu, loom.matrices){
  ## Ensure same cell order in UMAP embedding and in loom matrices
  seu.umap = Embeddings(seu, 'umap')
  geneset = rownames(seu)
  subset.cells <- intersect(rownames(seu.umap), colnames(loom.matrices$spliced))
  subset.genes <- intersect(geneset, rownames(loom.matrices$spliced))
  # subset and order
  loom.matrices = lapply(loom.matrices, function(loom) loom[subset.genes, subset.cells])
  if(length(subset.cells) < 2){ 
    loom.matrices = lapply(loom.matrices, function(loom) as(as.matrix(loom), "dgCMatrix"))
  }
  
  return(loom.matrices)
}


#4. Calculate cell distances matrix
calculate.cell.dist <- function(seu){
  seu.pca = Embeddings(seu, 'pca')
  cell.dist <- as.dist(1-armaCor(t(seu.pca), nthreads=8))
  cell.dist
}


# 5. Gene relative velocity
gene.velocity.estimates <- function(loom.matrices, cell.dist, fit.quantile){
  ## Compute gene relative velocity estimates
  velo.est <- gene.relative.velocity.estimates(loom.matrices$spliced, 
                                               loom.matrices$unspliced, 
                                               deltaT=1, 
                                               kCells=20, 
                                               cell.dist=cell.dist, 
                                               fit.quantile=fit.quantile, 
                                               n.cores=8)
  velo.est
}


# 6. Cell color vector
color.vector <- function(seu, cl_annot){
  ## Obtain vector for cell coloring
  suppressPackageStartupMessages(require(scales, lib.loc = "/usr/local/lib/R/site-library/scales"))
  
  labels <- seu@meta.data[[cl_annot]]
  ncol <- length(unique(labels))
  x <- as.factor(labels)
  levels(x) <- 1:length(levels(x))
  x <- as.numeric(x)
  colvec <- hue_pal()(ncol)[x]
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


### Execute

# 0. Params
# Annots
annotations <-  gsub(opt$annots, pattern = " ", replacement = "")
## Split by comma
annotations <- unlist(strsplit( annotations , ","))

# Seurat object/s
seu = readRDS(opt$input_seurat)


# subset it to cell type and disease stage if specified
if (opt$subset_celltype !='All') {
  cat(paste0('\nSUBSETTING TO KEEP ', opt$subset_celltype, '\n'))
  Idents(seu) <- seu@meta.data$cell.annotation
  seu <- subset(seu, idents=opt$subset_celltype)
}
if (opt$subset_disease_stage !='All') {
  cat(paste0('\nSUBSETTING TO KEEP ', opt$subset_disease_stage, '\n'))
  Idents(seu) <- seu@meta.data$Disease.status
  seu <- subset(seu, idents=opt$subset_disease_stage)
}

print('*****************')
print('after subsetting, the seurat object size is:')
print(seu)
print('*****************')


# Loom files 
cr_paths = read.csv(opt$cr_paths, header = T)

# 1. Read all loom files
loom.list <- list()
for (i in 1:nrow(cr_paths)){
  sample_name = as.character(cr_paths$sample[i])
  print(sample_name)
  cr.path = cr_paths$cr_path[i]
  print(cr.path)
  print('')
  loom.list[[sample_name]] = read.loom.file(cr.dir = cr.path)
}


seu_obj = seu
seu_name = opt$sample_name

# 2. Rename loom matrices
rename.loom.list <- list()
## Refactor
seu_obj@meta.data[["orig.ident"]] <- factor(seu_obj@meta.data[["orig.ident"]])

sample_name <- 'SLX-19940-SITTA2'
sample_name <- 'SLX-21151-SITTF7'
for(sample_name in levels(seu_obj@meta.data[["orig.ident"]])){
  print(sample_name)
  # seu object library suffix
  lib.suffixes = unlist(lapply(strsplit(colnames(seu_obj[, seu_obj@meta.data[["orig.ident"]] == sample_name]), "-"), function(x) x[2]))
  lib.suffix = paste0("-", levels(factor(lib.suffixes)))
  #print(lib.suffix)
  # rename loom matrix
  rename.loom.list[[sample_name]] = rename.col.outer(y = loom.list[[sample_name]], lib.suffix  = lib.suffix)
}

# Subset Seurat object to 20K cells if it is too big [to speed up, and to avoid crash when calculating cell-cell distances]
if(ncol(seu_obj) > 20000){
  if (!is.na(opt$barcode_file_path)) {
    barcodes_oI <- readRDS(opt$barcode_file_path)
    barcodes_oI <- barcodes_oI[barcodes_oI %in% rownames(seu_obj@meta.data)] 
    keep.barcodes <- c(barcodes_oI, sample(colnames(seu_obj), size = 20000-length(barcodes_oI)))
  } else {
      keep.barcodes <- sample(colnames(seu_obj), size = 20000)
  }
  seu_obj = seu_obj[, keep.barcodes]
}


# 3. Order and subset loom.matrices
loom.matrices.ord = lapply(rename.loom.list, function(loom) order.loom.cells(seu = seu_obj, loom.matrices = loom))
## 3.2. Filter matrices with < 2 cells
loom.matrices.ord = loom.matrices.ord[unlist(lapply(loom.matrices.ord, function(loom) ncol(loom$spliced) >= 2))]
## 3.3. Matrices must have the same number of rows
rownames.list = lapply(loom.matrices.ord, function(loom) rownames(loom$spliced))
intersect.rows = Reduce(intersect, rownames.list)
## 3.4. Subset
loom.matrices.ord = lapply(loom.matrices.ord, function(loom) lapply(loom, function(mat) mat[intersect.rows, ] ))


# 4. Concatenate loom matrices
## spliced
spliced.list = lapply(loom.matrices.ord, function(loom) loom$spliced)
spliced = Reduce(cbind, spliced.list)
## unspliced
unspliced.list = lapply(loom.matrices.ord, function(loom) loom$unspliced)
unspliced = Reduce(cbind, unspliced.list)

### 6. Subset seurat object
gene.subset = rownames(spliced)
cell.subset = colnames(spliced)
seu_obj <- seu_obj[ gene.subset, cell.subset]

# 7. Cell distances
cell.dist = calculate.cell.dist(seu = seu_obj)

# 8. Velocity estimates
velo.est = gene.velocity.estimates(loom.matrices = list("spliced" = spliced, "unspliced"=unspliced),
                                   cell.dist = cell.dist,
                                   fit.quantile = 0.2)
## save velocity estimates
saveRDS(velo.est,  paste0(opt$output_path, seu_name, "-velocitity.estimates.rds"))

# debugging
velo.est <- readRDS(paste0(opt$output_path, seu_name, "-velocitity.estimates.rds"))

# Iterate across annotations to use for coloring
for(cl_annot in annotations){
  print(cl_annot)
  
  # 6. Color vector
  col.vec = color.vector(seu = seu_obj, cl_annot = cl_annot)

  # 7. Velocity map on embedding
  print('plotting velocity map...')
  velocity.on.embedding.plot(seu = seu_obj, 
                             reduction='umap_harmony_t.0.1',
                             pdf_path = paste0(opt$output_path, seu_name, "-Velocity-Plot", cl_annot, ".pdf"),
                             grid_vals = c(20, 40, 60, 80), 
                             velo.est = velo.est, 
                             col_vec = col.vec,
                             sample_name = seu_name, cl_annot = cl_annot)
}



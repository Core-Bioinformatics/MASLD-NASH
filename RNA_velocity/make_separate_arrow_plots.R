rm(list=ls())

library(Seurat)
library(scales) # for hue_pal() in get colour vector
library(velocyto.R)
source('/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/7.RNA-Velocity/Scripts/make_separate_arrow_plots_functions.R')


# set params -------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
input_seurat_file <- args[1] # path to the seurat obj
reduction <- args[2] # the name of the reduction to use for plotting
input_chol_velObj_file <- args[3] # path to the velObjFile
input_hep_velObj_file <- args[4]
plotPath <- args[5] # path to output plot
filter.to.endstage <- args[6] # "filter" or anything else not to.

# run it -----------------------------------------------------------------------
print('***********************************************************************')
print(input_seurat_file)
print(reduction)
print(input_chol_velObj_file)
print(input_hep_velObj_file)
print(plotPath)
print(filter.to.endstage)
print('***********************************************************************')

print('loading Seu...')
seu <- readRDS(input_seurat_file)
print('loading Chol...')
velo.est.chol <- readRDS(input_chol_velObj_file)
print('loading Hep...')
velo.est.hep <- readRDS(input_hep_velObj_file)

if (filter.to.endstage=='filter') {
  print('filtering...')
  seu.filt <- subset(seu, subset=Disease.status=='end stage')
  seu <- seu.filt
}

# make cell color vector 
cl_annot='cell.annotation'
col_vec = color.vector(seu = seu, cl_annot = cl_annot )

# get cell position embeddings
emb = Embeddings(seu, reduction)

print('Making plot ....')
pdf(plotPath, width=10, height=10) 
  for (grid.val in c(10, 20, 30)) {
    make_embedding_plot(emb, cell.colors=col_vec, 
                        cell.border.alpha=0)  
    
    my.add.arrows(emb, 
                  velo.est.chol,
                  n=300,
                  scale="sqrt",
                  cell.colors = ac(col_vec, alpha=0.5),
                  arrow.color='#EC2055',
                  cex=0.8,
                  arrow.scale=5,
                  show.grid.flow=TRUE,
                  min.grid.cell.mass=0.5,
                  grid.n=grid.val, # determined number of arrows!
                  arrow.lwd=2,
                  do.par=F,
                  n.cores = 8,
                  cell.border.alpha = 0.1,
                  main='')
    
    my.add.arrows(emb, 
                  velo.est.hep,
                  n=300,
                  scale="sqrt",
                  cell.colors = ac(col_vec, alpha=0.5),
                  arrow.color='dodgerblue3',
                  cex=0.8,
                  arrow.scale=5,
                  show.grid.flow=TRUE,
                  min.grid.cell.mass=0.5,
                  grid.n=grid.val, # determined number of arrows!
                  arrow.lwd=2,
                  do.par=F,
                  n.cores = 8,
                  cell.border.alpha = 0.1,
                  main='')
  
  }
dev.off()

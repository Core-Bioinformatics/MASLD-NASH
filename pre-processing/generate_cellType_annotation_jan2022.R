rm(list=ls())
library(data.table)
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)

# Goal:
# use expression of marker genes per cluster to manually annotate the Jan2022
# dataset with celltypes


# Functions ----
mega_melt <- function(df, source) {
  df.m <- melt(df, measure.vars=names(df))
  names(df.m) = c('cellType', 'marker.gene')
  df.m$source = source
  
  return(df.m)
}


get_marker_expression_info <- function(all.markers, seu) {
  # for each cell barcode, calculate the number (and proportion) of the 
  # markers for each celltype expressed at a non-0 level in that cell. 
  
  count_expressed <- function(x) {
    return(sum(x > 0))
  }
  
  curr.cellType='Neutrophils'
  dt.list <- list()
  for (curr.cellType in unique(all.markers$cellType)) {
    print(curr.cellType)
    # get expression of the marker genes for the current cell type
    curr.markers <- all.markers[all.markers$cellType==curr.cellType,]
    curr.rows <- (rownames(seu@assays$RNA@counts) %in% curr.markers$marker.gene) # this has the full gene set
    curr.expression <- seu@assays$RNA@counts[curr.rows,]
    
    # if more than 1 marker gene for the cellType
    if (sum(curr.rows) > 1) {
      curr.num.markers.expressed <- apply(curr.expression, 2, count_expressed)
      # get the proportion of all the markers for the current cell type expressed in the cell
      proportion.of.markers.expressed <- curr.num.markers.expressed / nrow(curr.markers)
    } else { # if only 1 marker gene
      curr.num.markers.expressed <- curr.expression > 0
      proportion.of.markers.expressed <- curr.num.markers.expressed / 1
    }
    
    curr.dt <- data.table(data.frame('cell.id'= names(proportion.of.markers.expressed),
                                     'cellType' = curr.cellType,
                                     'expressed.marker.num'=curr.num.markers.expressed,
                                     'expressed.marker.proportion'=proportion.of.markers.expressed))
    dt.list[[curr.cellType]] <- curr.dt
  }
  marker.expr <- do.call('rbind', dt.list)
  return(marker.expr)
}

# cellannotation <- fixed_annotation[fixed_annotation$cell.id=='AAAGGTACATACCAGT-23',]
# collapse_cellTypes(cellannotation)
collapse_cellTypes <- function(cellannotation) {
  cellannotation <- cellannotation[order(cellannotation$cellType),]
  label <- paste0(cellannotation$cellType, collapse='//')
  return(label)
}

#annotation <- fixed_annotation
make_annotation_CT_collision_plot <- function(annotation, prop.th) {
  collisions <- annotation[, .('conflict'=collapse_cellTypes(.SD)), by=.(cell.id)]
  
  coll.counts <- collisions[, .('counts'=.N), by=.(conflict)]
  coll.counts <- coll.counts[order(coll.counts$counts, decreasing = T),]
  coll.counts$conflict <- factor(coll.counts$conflict, levels=unique(coll.counts$conflict))
  
  coll.counts <- coll.counts[coll.counts$counts > 5,] # don't sweat the small stuff!
  
  p <- ggplot(coll.counts, aes(x=conflict, y=counts))+
    geom_histogram(stat='identity')+
    scale_y_continuous(trans='log10')+
    xlab('annotation conflict')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=4))+
    ggtitle(paste0('proportion.th=', prop.th))
  return(p)
}

make_annotation_CT_distribution_plot <- function(annotation, prop.th) {
  n.cells <- length(unique(marker.expr$cell.id))
  n.annotated <- length(unique(annotation$cell.id))
  uniq.annotation <- annotation[annotation$num.putative.cellTypes==1,]
  n.uniq.annotated <- length(uniq.annotation$cell.id)
  
  type.counts <- uniq.annotation[, .('counts'=.N), by=.(cellType)]
  
  p <- ggplot(annotation, aes(x=cellType))+
    geom_histogram(stat='count', alpha=0.5)+
    geom_histogram(data=uniq.annotation, stat='count') +
    geom_text(data=type.counts, aes(x=cellType, y=counts, label=counts), vjust=1, color='white',
              size=1.5) +
    scale_y_continuous(trans='log10') +
    theme_bw()+
    xlab('')+
    ggtitle(paste0('proportion.th=', prop.th, '\n',
                   'n.annotated=', n.annotated, '\n',
                   'n.uniq.annotated=', n.uniq.annotated))+
    theme(axis.text = element_text(size=6),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          plot.title = element_text(size=6))
  
  return(p)
}

update_annotation <- function(annotation, meta, name, clusterNums, clusterCol) {
  update.rows <- meta[[clusterCol]] %in% clusterNums
  annotation[update.rows] <- name
  return(annotation)
}

filter_out_Wa_data <- function(seu) {
  seu@meta.data$keep <- !(seu@meta.data$orig.ident %in% c('Wa-1-postSort',
                                                          'Wa-1-preSort'))
  #m <- seu@meta.data
  #unique(m[, c('orig.ident', 'Patient.ID', 'manuscript.expt', 'keep')])
  seu.filt <- subset(seu, subset=keep==TRUE)
  seu.filt@meta.data$keep <- NULL
  #m <- seu.filt@meta.data
  #unique(m[, c('orig.ident', 'Patient.ID', 'manuscript.expt')])
  return(seu.filt)
}

seurat_preprocess <- function(seu){

  # Normalizarion
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


# MAIN ----

# path to preprocessed seurat object, (produced by do_preprocessing.sh)
seuObjDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/6.Aggr-Jan2022/'
seuObjPath <- paste0(seuObjDir, 'Aggr_Jan2022_pre-processed.rds')
outObjFile <- 'Aggr_Jan2022_annotated.rds'
seuObjOutPath <- paste0(seuObjDir, outObjFile)

# path to plot output directory
plotDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/6.Aggr-Jan2022/cellType_annotation/'
dir.create(plotDir)

#1. load the Seurat object
seu <- readRDS(seuObjPath)
seu

## 2. Load Marker Genes
markerGeneDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/2.Markers-expression/'

# now ONLY use the most discriminative markers Ruben previously identified, 
# and the B-cell ones
all_the_markers <- fread(paste0(markerGeneDir, "0.Marker_genes_sets/marker_genes_v4.tsv"), sep = "\t")
am_m <- mega_melt(all_the_markers, 'AM')
am_m$cellType <- as.character(am_m$cellType)
am_m$cellType[am_m$cellType=='Stellate cells'] <- 'Stellate'

all.markers <- am_m
rm(all_the_markers, am_m)


# 3. tidy up all.markers:
# remove duplicates
all.markers <- all.markers[!(duplicated(all.markers[, c('cellType', 'marker.gene')])), ]
# remove blank genes (from the csv format)
all.markers <- all.markers[all.markers$marker.gene != '']

# remove any putative markers not detected in the data

### plot expression of marker genes by the clusters of cells identified in pre-processing
resns = c(2.0)
r = 2.0
for (r in resns) {
  clustercol = paste0('SCT_snn_res.', r)
  
  print(clustercol)
  p.umap <- DimPlot(seu, reduction='umap', group.by=clustercol)
  p.umap.split <- DimPlot(seu, reduction='umap', 
                          group.by=clustercol,split.by = clustercol, ncol=7)
  p.heatmap <- DoHeatmap(seu, features=all.markers$marker.gene, group.by=clustercol)
  
  ggsave(paste0(plotDir,
                'marker_umap_', clustercol, '.pdf'), 
         plot = p.umap, 
         width=20, height=20)
  ggsave(paste0(plotDir, 
                'marker_umap_split_', clustercol, '.pdf'),
         plot=p.umap.split,
         width=20, height=20)
  ggsave(paste0(plotDir, 'marker_heatmap_', clustercol, '.pdf'),
         plot=p.heatmap,
         width=30, height=30)
}

# based on the above plots, assign the clusters to cellTypes ----
# A FIRST PASS TO ASSIGN CELLS BASED ON CLUSTERS
# based on marker by cluster in "marker_heatmap_SCT_snn_res.2_annotated.pdf"
hep.clusters = c(0,1,2,3,4,5,6,7,8,10,11,12,17,18,19,22, 27, 29, 31, 34, 39, 43) # 19 half hep, halve Lymph - will need refining
chol.clusters = c(15, 24, 28, 42)
stellate.clusters = c(23, 30,33)
endothelial.clusters = c(9,16,20,21, 32, 40, 41)
lympho.clusters = c(13,14, 36, 38)
macro.clusters = c(25, 26, 37) # half of 25 is neutrophils 
Bcell.clusters = c(35)
neutro.clusters = c()
unassigned.clusters = c()

# assign cellType annotations to seurat object
meta <- seu@meta.data
annotation = rep('notDone', nrow(meta))
annotation <- update_annotation(annotation, meta, 'Hepatocytes', hep.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Cholangiocytes', chol.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Stellate', stellate.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Endothelial', endothelial.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Lymphocytes', lympho.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Macrophages', macro.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'B-cell', Bcell.clusters, 'SCT_snn_res.2')
annotation <- update_annotation(annotation, meta, 'Neutrophil', neutro.clusters, 'SCT_snn_res.2')
# unannotated labelled with clusters
annotation[annotation=='notDone'] <- as.character(seu@meta.data$SCT_snn_res.2[annotation=='notDone'])
#annotation <- update_annotation(annotation, meta, 'unknown', unassigned.clusters, 'SCT_snn_res.2')

# add the annotation to the seu obj
seu <- AddMetaData(seu, 'Annot.Jan2021.tmp', col.name='cell.annotation.version')
seu <- AddMetaData(seu, annotation, col.name='cell.annotation')

unique(seu@meta.data[, c('cell.annotation', 'SCT_snn_res.2')])

# sanity plot how it looks currently
# sanity plot UMAP 
p <- DimPlot(seu, reduction='umap', group.by='cell.annotation')
psplit <- DimPlot(seu, reduction='umap', group.by='cell.annotation', split.by='cell.annotation', ncol=3)
p.both <- DimPlot(seu, reduction='umap', group.by='SCT_snn_res.2', split.by='cell.annotation', ncol=3)

ggsave(paste0(plotDir, 'cellType_UMAP.pdf'), 
       plot=p, 
       width=10, height=10)
ggsave(paste0(plotDir, 'cellType_UMAP_facetted.pdf'), 
       plot=psplit, 
       width=20, height=20)

p.list <- list()
for (ct in unique(seu@meta.data$cell.annotation)) {
  Idents(seu) <- seu@meta.data$cell.annotation
  seu.sub <- subset(seu, idents=ct)
  p <- DimPlot(seu.sub, reduction='umap', group.by='SCT_snn_res.2')+
    ggtitle(ct)
  p.list[[ct]] <- p
}
p.all <- plot_grid(plotlist=p.list)
ggsave(paste0(plotDir, 'cellType_UMAPgroup&celltype.pdf'), 
       plot=p.all, 
       width=20, height=20)

# tidy the clusters that look like have multiple celltypes:
meta <- seu@meta.data
umap.embeddings <- seu@reductions$umap@cell.embeddings
expr <- seu@assays$SCT@scale.data

# fixing lymphocytes
meta$cell.annotation[meta$cell.annotation=='Lymphocytes' &
                       umap.embeddings[, 1] < -5 &
                       umap.embeddings[, 2] < 5] <- 'Hepatocytes'
# B-cell
# label B2 popn
meta$cell.annotation[meta$cell.annotation=='B-cell' &
                       (umap.embeddings[, 2] < 1 &
                          umap.embeddings[, 1] > -5)] <- 'B-cell 2'
# label little ones in between hep and chol
meta$cell.annotation[meta$cell.annotation=='B-cell' &
                       (umap.embeddings[, 2] < 3 &
                          umap.embeddings[, 1] < -6)] <- 'Hepatocytes'
# label B1 popn
meta$cell.annotation[meta$cell.annotation=='B-cell' &
                       (umap.embeddings[, 2] > 5 & umap.embeddings[, 1]< 3)] <- 'B-cell 1'
# label the ones going up to B1
meta$cell.annotation[meta$cell.annotation=='B-cell 1' & 
                       (umap.embeddings[, 1] < -2 &
                          umap.embeddings[, 2] < 11)] <- 'unknown'
meta$cell.annotation[meta$cell.annotation=='B-cell'] <- 'unknown'

# fixing Hepatocytes 
# which are really lymphocytes
meta$cell.annotation[meta$cell.annotation=='Hepatocytes' &
                       (umap.embeddings[, 1] < 5 &
                          umap.embeddings[, 1] > -4 &
                          umap.embeddings[, 2] > 2.5 &
                          umap.embeddings[, 2] < 7.5)] <- 'Lymphocytes'

# which are in hep UMAP group, but lymph ofr stellate clusters
meta$cell.annotation[meta$cell.annotation %in% c('Stellate', 'Lymphocytes', 'Neutrophils') &
                       (umap.embeddings[, 1] > 0 &
                          umap.embeddings[, 2] > -7.5 &
                          umap.embeddings[, 2] < 1)] <- 'Hepatocytes'

# fixing Macrophages
meta$cell.annotation[meta$cell.annotation=='Macrophages' &
                       expr['FCN1',] > 0] <- 'Neutrophils'

seu <- AddMetaData(seu, meta$cell.annotation, col.name='cell.annotation.refined')

# sanity plot how it looks currently
# sanity plot UMAP
p <- DimPlot(seu, reduction='umap', group.by='cell.annotation.refined')
psplit <- DimPlot(seu, reduction='umap', group.by='cell.annotation.refined', split.by='cell.annotation.refined', ncol=3)
p.both <- DimPlot(seu, reduction='umap', group.by='SCT_snn_res.2', split.by='cell.annotation.refined', ncol=3)

ggsave(paste0(plotDir, 'cellType_UMAP_refined.pdf'),
       plot=p,
       width=10, height=10)
ggsave(paste0(plotDir, 'cellType_UMAP_facetted_refined.pdf'),
       plot=psplit,
       width=20, height=20)

# check how many cells of each type identified
table(seu@meta.data$cell.annotation.refined)
 
# only keep the refined annotation
seu@meta.data$cell.annotation <- seu@meta.data$cell.annotation.refined
seu@meta.data$cell.annotation.refined <- NULL

# sanity plot UMAP 
p <- DimPlot(seu, reduction='umap', group.by='cell.annotation')
psplit <- DimPlot(seu, reduction='umap', group.by='cell.annotation', split.by='cell.annotation', ncol=3)

# 5. save the cellType annotated Seuarat object
saveRDS(seu,
        file=seuObjOutPath)

# save split by disease state for use in cellPhoneDB
#seu <- readRDS(file=seuObjOutPath)
m <- seu@meta.data
dis.list <- SplitObject(seu, split.by='Disease.status')
dis <- names(dis.list)[1]
for (dis in names(dis.list)) {
  print(dis)
  curr.seu <- dis.list[[dis]]
  
  dis.safe <- str_replace_all(dis, ' ', '-')
  dis.safe <- str_replace_all(dis.safe, '/', '')
  print(dis.safe)
  out.path <- paste0(seuObjDir, 'Aggr_by_disease/', dis.safe, '_', outObjFile)
  saveRDS(curr.seu, file=out.path)
}


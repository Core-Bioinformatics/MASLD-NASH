rm(list=ls())

library(Seurat)
library(stringr)


# reads clustering info from seurat object "cluster_seu_name" and adds it to 
# "projected_seu"

projDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'
seuDir <- '6.Aggr-Jan2022/'

cluster_seu_name <- 'Aggr_Jan2022_Chol_harmony_th=0.rds' # the seurat info with the cluster is interested in
projected_seu_name <- 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds' # the seurat object going to use for the UMAP

cluster_seu <- readRDS(paste0(projDir, seuDir, cluster_seu_name))
projected_seu <- readRDS(paste0(projDir, seuDir, projected_seu_name))

cluster_col <- 'SCT_snn_harmony_t.0.0.4'
cluster_oI <- '9'
new_cluster_name <- 'chol-9'

# add cluster info to cell type column of projected_seu
barcodes <- rownames(cluster_seu@meta.data)[cluster_seu@meta.data[[cluster_col]]==cluster_oI]

projected_seu@meta.data$cell.annotation[rownames(projected_seu@meta.data) %in% barcodes] <- new_cluster_name
unique(projected_seu@meta.data$cell.annotation)


outFileName <- str_replace_all(projected_seu_name, '.rds', '')
saveRDS(projected_seu, 
        file=paste0(projDir, seuDir, 'cluster_of_interest_objects/',
              outFileName, '_', new_cluster_name, '.rds'))

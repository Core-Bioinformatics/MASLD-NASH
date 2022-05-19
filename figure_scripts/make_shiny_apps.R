rm(list=ls())
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ShinyCell))
library(stringr)
library(readr)


# Control where input is ----
seu_dir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/6.Aggr-Jan2022/'    

# Control where output goes ----
outDir <- '/home/USSR/awc30/liver_project/manuscript_figures/shiny_public/'
dir.create(outDir, recursive=T)

cellTypes <- list('B-cell1'='B-cell.1', 
                  'B-cell2'='B-cell.2', 
                  'Chol'='Cholangiocytes', 
                  'Endo'='Endothelial', 
                  'Lymph'='Lymphocytes', 
                  'Macro'='Macrophages', 
                  'Neutro'='Neutrophils',
                  'Stell'='Stellate', 
                  'Hep'='Hepatocytes')




celltype <- names(cellTypes)[3]
for (celltype in names(cellTypes)) {
  print(celltype)
  
  # get and tidy input ---
  seuFile <- paste0('Aggr_Jan2022_', celltype, '_harmony_th=0.rds')
  curr.seu <- readRDS(paste0(seu_dir, seuFile))
  m <- curr.seu@meta.data
  
  # sanitize the meta-data
  del.cols <- c('manuscript.expt',
                'Steatosis', 'Ballooning', 'Inflammation', 'Fibrosis.score..F0.4.',
                'BMI', 'Ethnic.group', 'Diabetes.type.2', 'hypertension', 
                'dyslipidaemia', 'cardiovascular.disease', 'obstructive.sleep.apnoea',
                'percent.mt.RNA', 'percent.rp.RNA', 'nCount_SCT', 'nFeature_SCT', 
                'seurat_clusters', 'cell.annotation.version')
  del.clust.cols <- c(grep('SCT_snn_res', names(curr.seu@meta.data), value=T),
                      grep('SCT_snn_orig', names(curr.seu@meta.data), value=T))
  del.cols <- c(del.cols, del.clust.cols)
  
  
  for (del.col in del.cols) {
    curr.seu@meta.data[[del.col]] <- NULL
  }

  # make the app
  f_out <- cellTypes[[celltype]]
  scConf <- createConfig(curr.seu)
  makeShinyApp(curr.seu,
               gex.assay='SCT', 
               gex.slot='data',
               scConf,
               gene.mapping=T,
               shiny.dir=paste0(outDir, f_out),
               shiny.title=f_out)
  
}




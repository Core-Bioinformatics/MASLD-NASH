rm(list=ls())

library(Seurat)
library(harmony)
library(mcclust)
library(cowplot)
library(ggplot2)
library(data.table)
library(stringr)

myDimPlot <- function(meta.data, group.by) {
  p <- ggplot(meta.data, aes_string(x='UMAP1', y='UMAP2', color=group.by))+
    geom_point(size=0.01)+
    theme_classic()
  return(p)
}


seurat_preprocess <- function(seu){
  # Normalizarion
  seu <- SCTransform(seu, assay = "RNA", variable.features.n = 3000) # variable.features.n set to default
  
  # Scale Data (including all features)
  seu <- ScaleData(seu, assay = "SCT", features=NULL, do.scale = TRUE, do.center = TRUE)
  
  # Run PCA
  seu <- RunPCA(seu, assay = "SCT", verbose = F)
  return(seu)
}


do_orig_UMAP_and_clustering <- function(seu, cluster.resns) {
  # original as in using PCA, rather than harmony corrected PCA
  
  seu <- FindNeighbors(seu, 
                       reduction = 'pca') # neighbours are used in FindClusters, but not in RunUmap
  
  seu <- RunUMAP(seu, 
                 reduction='pca',
                 dims = 1:50, 
                 reduction.name='umap')
  
  for (c.res in cluster.resns) {
    seu <- FindClusters(seu, 
                        resolution = c.res,
                        alg=1)
    
    seu <- AddMetaData(seu, 
                       seu@meta.data$seurat_clusters, 
                       col.name = paste0('SCT_snn_', 'orig', '.', c.res))
  }
  return(seu)
}


harmony_and_reprocess <- function(seu, theta, reduction.name, cluster.resns) {
  # theta is theta for harmony (strength of batch correction)
  # reduction.name is name of reduction to give to harmony reduction (and will 
  # be used in naming UMAP reduction, and clusters)
  # cluster.resns is vector of clustering resolutions to cluster on
  
  # Run Harmony --- 
  seu <- RunHarmony(seu, 
                    group.by.vars = 'Patient.ID',
                    reduction='pca',
                    theta=theta, 
                    reduction.save=reduction.name)
  
  seu <- FindNeighbors(seu, 
                       reduction = reduction.name)
  
  
  # UMAP - runs on the reduced dims directly - not on neighborhood graph
  seu <- RunUMAP(seu, 
                 reduction=reduction.name,
                 dims = 1:50, 
                 reduction.name=paste0('umap_', reduction.name),
                 )
  
  # Clusters
  for (c.res in cluster.resns) {
    # uses the Neighbors last found
    seu <- FindClusters(seu, 
                        resolution = c.res, 
                        alg = 1
                        )
    seu <- AddMetaData(seu, 
                       seu@meta.data$seurat_clusters, 
                       col.name = paste0('SCT_snn_', reduction.name, '.', c.res))
  }
  
  return(seu)
}


plot_cluster_composition <- function(m, by.col, clust.col, titleStr, show.legend) {
  
  # calculate proportions of each cluster is each group
  m <- data.table::data.table(m)
  plt.dt <- m[, .(counts=.N), by=.(get(by.col), get(clust.col))]
  names(plt.dt)[names(plt.dt)=='get'] <- by.col
  names(plt.dt)[names(plt.dt)=='get.1'] <- clust.col
  
  plt.dt[, proportion:=counts / sum(counts), by=.(get(clust.col))]
  
  # make barplots coloured by by.col with shared legend
  p1 <- ggplot(plt.dt, aes_string(x=clust.col, y='counts', fill=by.col))+
    geom_bar(stat='identity')+
    theme_bw()+
    theme(legend.position='top')
  
  p2 <- ggplot(plt.dt, aes_string(x=clust.col, y='proportion', fill=by.col))+
    geom_bar(stat='identity')+
    theme_bw()+
    theme(legend.position='none')
  
  # get shared legend
  legend <- get_legend(p1)
  p1 <- p1+theme(legend.position='none')
  
  p.bars <- plot_grid(p1, p2, nrow=1)
  if (show.legend) {
    p.leg <- plot_grid(legend, p.bars, ncol=1, rel_heights=c(0.2, 1))
  } else {
    p.leg <- p.bars
  }
  
  # plot entropy per cluster
  entropy.dt <- plt.dt[, .(entropy=-1*sum(proportion * log2(proportion))), by=.(get(clust.col))]
  names(entropy.dt)[names(entropy.dt)=='get'] <- clust.col
  
  p3 <- ggplot(entropy.dt, aes_string(x=clust.col, y='entropy'))+
    geom_bar(stat='identity', position=position_dodge2())+
    theme_bw()
  
  p.all <- plot_grid(p.leg, p3, nrow=1, rel_widths = c(1, 0.3))
  
  title <- ggdraw()+draw_label(titleStr)
  p.finished <- plot_grid(title, p.all, ncol=1, rel_heights = c(0.1, 1))
  
  return(p.finished)
}


plot_orig_and_harm_cluster_purities <- function(m, by.col, c.resns, show.legend) {
  plot.list <- list()
  # c.c=1.6
  for (c.c in c.resns) {
    c.cols <- names(m)[grepl(paste0(c.c, '$'), names(m)) & 
                         !(grepl('_res', names(m)))]
    if (length(c.cols) != 2) {
      print('columns:')
      print(c.cols)
      stop('Regex didnt find the expected 2 columns!')
    }
    orig.col <- c.cols[grepl('_orig', c.cols)]
    harm.col <- c.cols[grepl('harmony', c.cols)]
    
    p.orig <- plot_cluster_composition(m, by.col, orig.col, paste0('Orig.clusts. c.res=', c.c), show.legend)
    p.harm <- plot_cluster_composition(m, by.col, harm.col, paste0('Harmony.clusts. c.res=', c.c), show.legend)
    plot.list[[paste0(c.c)]] <- plot_grid(p.orig, p.harm, nrow=1)
  }
  p.grid <- plot_grid(plotlist=plot.list, ncol=1)
  return(p.grid)
}

#str <- 'NASH w/o cirrhosis'
make_names_safe_to_plot <- function(str) {
  out <- str
  out <- str_replace_all(out, '/', '')
  out <- str_replace_all(out, ' ', '_')
}


# goal - to apply harmony to get patients within a dataset lining up, without
# making all the disease states collapsing



#dataset.to.use = args[1]


### INPUT FILES-----
projDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'

dataset.to.use = 'CG_jan2022'

if (dataset.to.use == 'CG_only') {
  print('CG only data')
  seuDir <- '6.Aggr-Dec2021_v2_CG_only/' # which data version to use
  seuFile <- 'Aggr_Dec2021_CGonly_annotated.rds'
} else if (dataset.to.use == 'CG+Wa') {
  print('CG + Wa data')
  seuDir <- '6.Aggr-Dec2021_v2_CG+Wa/'
  seuFile <- 'Aggr_Dec2021_CG+Wa_annotated.rds'
} else if (dataset.to.use == 'CG+id98') {
  seuDir <- '6.Aggr-Dec2021_v2_CG+Wa+id98/'
  seuFile <- 'Aggr_Dec2021_CG+id98_annotated.rds'
} else if (dataset.to.use == 'CG_jan2022') {
  seuDir <- '6.Aggr-Jan2022/'
  seuFile <- 'Aggr_Jan2022_annotated.rds'
} else if (dataset.to.use == 'CG_jan2022_filt') {
  seuDir <- '6.Aggr-Jan2022_filt_prolif/'
  seuFile <- 'Aggr_Jan2022_no_prolif_annotated.rds'
} else {
  stop('invalid dataset.to.use')
}


### RUN PARAMETERS ------
args = commandArgs(trailingOnly = T)
cellStr <- args[1] # 'All', 'Chol', 'Hep' 'Chol-Hep'
sep_by_dis <- args[2] # "T" or "F"
theta_in <- args[3]

#all.thetas <- '2' # args # 2 # theta parameter to harmony (regularisation strength)

# which version of the harmony correction analysis to run 
# v1 (don't reprocess on subsetting), 
# v1.5 (reprocess when subset to cells) 
# or v2 (reprocess when subset to cells and when subset to disease state)
version <- 'v1.5' 
#all.thetas <- c('0', '0.1', '1', '2', '4')
all.thetas <- c(theta_in)
all.thetas <- sapply(all.thetas, as.numeric)
clust.resns <- c(0.1, 0.4, 0.8, 1.2, 1.6)

# if apply harmony seperately to each disease state, then recombine
# seperate_by_disease <- T # FALSE or TRUE
if (sep_by_dis=='T') {
  print('sep_by_dis==T!!!')
  seperate_by_disease <- T # FALSE or TRUE
} else if (sep_by_dis=='F') {
  seperate_by_disease <- F
} else {
  stop('inadequate sep_by_dis specified.')
}
# if seperate_by_disease==T, specify the disease stages to analyse if don't want to run for all
# NULL will run for all
dis.stages.to.analyse <- NULL # c('control') # c("NASH_wo_cirrhosis", "NAFLD", "Healthy_control", "Disease_control", "end_stage", "NASH_with_cirrhosis", "control")

# if subset to only particular celltypes
# args = commandArgs(trailingOnly = T)
# cellStr <- args[1] # 'Chol', or 'Hep', or 'Chol-Hep'
#cellStr <- 'Chol'


### OUTPUT FILES ------
# will be made as subdirs within seuDir
outDirectory <- paste0('harmony_', cellStr, '_', version, '_', dataset.to.use, '/') # where all the output of the current run will go
plotDirectory <- paste0(outDirectory, '/harmony_plots_', cellStr, '_only_sep_disease=', sep_by_dis, '/') # plots
clustDirectory <- paste0(outDirectory, '/harmony_clusters_', cellStr, '_only_sep_disease=', sep_by_dis, '/') # csvs of clusters
outSeuDirectory <- paste0(outDirectory, '/seuObjs_', cellStr, '_only_sep_disease=', sep_by_dis, '/')



# Shouldn't need to change anything below here -----

if (cellStr == 'Chol') {
  subsetCells <- c('Cholangiocytes') # c('Hepatocytes')  # NULL
} else if (cellStr =='Hep') {
  subsetCells <- c('Hepatocytes')
} else if (cellStr == 'Chol-Hep') {
  subsetCells <- c('Cholangiocytes', 'Hepatocytes')
} else if (cellStr == 'All'){
  subsetCells <- NULL
} else if (cellStr == 'Endo-Hep') {
  subsetCells <- c('Endothelial', 'Hepatocytes')
} else {
  stop('incorrect cellStr')
}

print('params:')
print(paste0('infile:\n', seuDir, seuFile))
print(paste0('out:\n', outDirectory))
print(paste0('sep.disease:\n', seperate_by_disease))
print(paste0('subset cells:\n', subsetCells))
print(paste0('version:\n', version))
print(paste0('thetas:\n', all.thetas))
print(paste0('clust.resns:\n', clust.resns))




theta <- all.thetas[2]
for (theta in all.thetas) {
  
  # name of the output seurat file
  outFile <- paste0('Aggr_Jan2022_filt_', cellStr, '_harmony_th=', theta, '.rds')
  
  
  print(paste0('running harmony with theta = ', theta))
  
  # setup paths---
  # input
  seuPath <- paste0(projDir, seuDir, seuFile)
  # outputs
  plotDir <- paste0(projDir, seuDir, plotDirectory)
  clustDir <- paste0(projDir, seuDir, clustDirectory)
  outDir <- paste0(projDir, seuDir, outSeuDirectory) # where the seurat objects saved
  
  dir.create(plotDir, recursive = T)
  dir.create(clustDir, recursive = T)
  dir.create(outDir, recursive=T)
  
  # load seurat data
  print(seuPath)
  seu <- readRDS(seuPath)
  
  unique(seu@meta.data[, c('Patient.ID', 'Disease.status')])
  
  # convert patient 98 from "Disease control" to "Healthy control" (is diseased), 
  # as Vas reckons looks ok based on Histology
  seu@meta.data$Disease.status[seu@meta.data$Patient.ID==98] <- 'Healthy control'
  
  #### subset to cell types if required ---
  Idents(seu) <- 'cell.annotation'
  if (!(is.null(subsetCells))){
    seu <- subset(seu, idents=subsetCells)
    
    # if subsetting cells, should probably re-preprocess (SCTTansform, variable genes, PCA)
    if (version %in% c('v1.5', 'v2')) {
      seu <- seurat_preprocess(seu)
    }
  }
  
  m <- seu@meta.data
  
  if (seperate_by_disease) {
    seu_list <- SplitObject(seu, split.by='Disease.status')
    names(seu_list) <- sapply(names(seu_list), make_names_safe_to_plot)
    
    # subset to only the diseae stages want to run for
    if (!(is.null(dis.stages.to.analyse))) {
      seu_list <- seu_list[dis.stages.to.analyse]
    }
    
    # repreprocess (SCTtransform, variable genes, PCA) within each disease stage
    if (version=='v2') {
      seu_list <- sapply(seu_list, seurat_preprocess)
      gc()
    }
  }
  
  #### run harmony & reprocess it
  reduction.name <- paste0('harmony_t.', theta)
  
  if (seperate_by_disease==F) {
    # get clusters using pcas and variable genes recalculated for the current 
    # cells
    seu <- do_orig_UMAP_and_clustering(seu, 
                                       clust.resns)
    
    # apply harmony, and recalculate the clusters
    seu <- harmony_and_reprocess(seu, theta, reduction.name,
                                 clust.resns)
  } else {
    for (n in names(seu_list)) {
      print(n)
      seu_list[[n]] <- do_orig_UMAP_and_clustering(seu_list[[n]], 
                                                   clust.resns)
      
      seu_list[[n]] <- harmony_and_reprocess(seu_list[[n]], theta, reduction.name,
                                             clust.resns)
      
    }
  }
  
  # summary plots of the harmony reduction
  if (seperate_by_disease==F) {
    p1 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='cell.annotation')
    p2 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Disease.status')
    p3 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Patient.ID')
    p4 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='manuscript.expt')
    clust.to.plot <- paste0('SCT_snn_harmony_t.', theta, '.0.4') # plot clustering for current theta, c.res=0.4
    p5 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by=clust.to.plot)
    p.all <- plot_grid(p1, p2, p3, p4, p5, ncol=2)
    ggsave(file=paste0(plotDir, 'summary_UMAPs_theta=', theta, '.png'),
           width=20, height=10)
    
    p <- DimPlot(seu, reduction=paste0('umap_', reduction.name), 
                  group.by='Disease.status', 
                  split.by='Patient.ID', 
                 ncol=7)
    ggsave(file=paste0(plotDir, 'summary_UMAPs_theta=', theta, '_split_patient.png'),
           width=10, height=10)
    
    # plot the mixing of patients per cluster for the current theta value
    m <- seu@meta.data
    c.resns <- clust.resns
    p.dis <- plot_orig_and_harm_cluster_purities(m, 'Disease.status', c.resns, show.legend=T)
    p.patient <- plot_orig_and_harm_cluster_purities(m, 'Patient.ID', c.resns, show.legend=F)
    
    ggsave(file=paste0(plotDir, 'cluster_purity_theta=', theta, '_disease.status.pdf'),
           plot=p.dis,
           width=15, height=20)
  
    ggsave(file=paste0(plotDir, 'cluster_purity_theta=', theta, '_patientID.pdf'),
           plot=p.patient,
           width=15, height=20)
    
    # sanity_umap
    p <- DimPlot(seu, reduction='umap', group.by=paste0('SCT_snn_harmony_t.', theta, '.0.1'))
    p.h <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by=paste0('SCT_snn_harmony_t.', theta, '.0.1'))
    p.grid <- plot_grid(p,p.h, nrow=1)
    ggsave(file=paste0(plotDir, 'cluster_purity_theta=', theta, '_sanity_UMAP.pdf'),)
    
    # write the clusters to csv, so can compare cluster similarities across theta thresholds
    # m <- seu@meta.data
    # m$barcodes <- rownames(m)
    # fwrite(m,
    #        paste0(clustDir, 'harmony_clusters_theta=', theta, '.csv'))
    gc()
    print('about to save SeuObj:')
    print(seu)
    
    saveRDS(seu, 
            file=paste0(outDir, outFile))
    
  } else {
    for (n in names(seu_list)) {
      seu <- seu_list[[n]]
      
      p1 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='cell.annotation')
      p2 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Disease.status')
      p3 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='Patient.ID')
      p4 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by='manuscript.expt')
      clust.to.plot <- paste0('SCT_snn_harmony_t.', theta, '.0.4') # plot clustering for current theta, c.res=0.4
      p5 <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by=clust.to.plot)

      p.all <- plot_grid(p1, p2, p3, p4, p5, ncol=2)
      ggsave(file=paste0(plotDir, 'summary_UMAPs_theta=', theta, '_', n, '.png'),
             width=20, height=10)
      
      p <- DimPlot(seu, reduction=paste0('umap_', reduction.name), 
                   group.by='Disease.status', 
                   split.by='Patient.ID', 
                   ncol=7)
      ggsave(file=paste0(plotDir, 'summary_UMAPs_theta=', theta, '_', n, '_split_patient.png'),
             width=10, height=10)
      
      # plot the mixing of patients per cluster for the current theta value
      m <- seu@meta.data
      c.resns <- clust.resns
      p.dis <- plot_orig_and_harm_cluster_purities(m, 'Disease.status', c.resns, show.legend=T)
      p.patient <- plot_orig_and_harm_cluster_purities(m, 'Patient.ID', c.resns, show.legend=F)
      
      ggsave(file=paste0(plotDir, 'cluster_purity_theta=', theta, '_', n, '_disease.status.pdf'),
             plot=p.dis,
             width=15, height=20)
      
      ggsave(file=paste0(plotDir, 'cluster_purity_theta=', theta, '_', n, '_patientID.pdf'),
             plot=p.patient,
             width=15, height=20)
      
      # sanity_umap
      p <- DimPlot(seu, reduction='umap', group.by=paste0('SCT_snn_harmony_t.', theta, '.0.1'))
      p.h <- DimPlot(seu, reduction=paste0('umap_', reduction.name), group.by=paste0('SCT_snn_harmony_t.', theta, '.0.1'))
      p.grid <- plot_grid(p,p.h, nrow=1)
      ggsave(file=paste0(plotDir, 'cluster_purity_theta=', theta, '_sanity_', n, '_UMAP.pdf'),
             width=20, height=10)
      
      
      # write the clusters to csv, so can compare cluster similarities across theta thresholds
      # m <- seu@meta.data
      # m$barcodes <- rownames(m)
      # data.table::fwrite(m,
      #        paste0(clustDir, 'harmony_clusters_theta=', theta, '_', n, '.csv'))
      
      # write seurat to file
      gc()
      print('about to save seu obj:')
      print(seu)
      saveRDS(seu, 
              file=paste0(outDir, outFile, '_', n, '.rds'))
      
    }
  }
}

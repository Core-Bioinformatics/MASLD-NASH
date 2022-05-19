rm(list=ls())
library(Seurat)
library(ggplot2)
library(data.table)
source('/home/USSR/awc30/liver_project/manuscript_figures/manuscript_UMAPs&dotplots/plot_functions.R')


projDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'
seuDir <- '6.Aggr-Jan2022/'
seuFile <- 'Aggr_Jan2022_Hep_harmony_th=0.rds'

plotdir <- '/home/USSR/awc30/liver_project/manuscript_figures/manuscript_UMAPs&dotplots/plots/'



seu <- readRDS(paste0(projDir, seuDir, seuFile))

disease.order <- c('Healthy control', 'NAFLD', 'NASH w/o cirrhosis',
                   'NASH with cirrhosis', 'end stage')
seu@meta.data$Disease.status <- factor(seu@meta.data$Disease.status,
                                      levels=disease.order)

#seu@reductions$

# 2A --------
L <- make_UMAP_plot(seu, 
                    'umap_harmony_t.0', 
                    'Disease.status')

L[[1]]
pdf(file=paste0(plotdir, 'fig2A.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()


# 2B : DEG bubble plot and 3 DEG in each disease status

# manually get most signif DE genes previously calcualted from, which are 
# specific to a stage
# '/home/USSR/awc30/liver_project/DEG_analysis_disease/results/Dis.vs.all/DE_results/Hepatocytes/'
healthy.up <- c('CRP', 'EVL', 'FGA') #c('CRP', 'GC', 'FGA')
nafld.up <- c('CYP2C19', 'CTNNA3', 'ALDH1L1') #c('HMGCS2', 'RORC', 'ALDH1L1')
nash.up <- c('ZNF385D', 'AP002026.1', 'MLIP')#c('CFAP57',
          #   'PRKCE',
          #   'ZNF385D')
nash.cirr.up <- c('LINC00486', 'HNF1A-AS1', 'FTCD') #c('LINC00486','NECAB2','APOC1')
end.up <-c('CYP3A7', 'SERPINE1', 'IGFBP1') #c('ERRFI1', 'PLPP3', 'JUN')

plt.genes <- c(healthy.up, nafld.up, nash.up, nash.cirr.up, end.up)
Idents(seu) <- seu@meta.data$Disease.status
p <- DotPlot(seu, assay='SCT',
             features=plt.genes, 
             scale=T)
p <- p+coord_flip()+
  xlab('')+
  ylab('')+
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p
p.nolab <- p + theme(legend.position = 'none')

pdf(file=paste0(plotdir, 'fig2B.pdf'), height=7, width=5)
  print(p)
  print(p.nolab)
dev.off()



# 2C : UMAPs for expression of the corresponding genes
plist <- list()
mg <- plt.genes[3]
for (mg in plt.genes) {
  p <- FeaturePlot(seu, reduction='umap_harmony_t.0', features=mg, slot='data')
  # p <- p+scale_color_viridis_c()+
  p <- p+ xlab('UMAP1') +
    ylab('UMAP2')+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank())
  plist[[mg]] <- p
}

pdf(file=paste0(plotdir, '2C_DEG_UMAPs.pdf'), height=7, width=7)
for (mg in names(plist)) {
  print(plist[[mg]])
}
dev.off()




# Do the same UMAPs and DEGs across disease for each of the other cell types
# for Supplemental figure 2 ----------------------------------------------------

# DE results already computed to:
# '/home/USSR/awc30/liver_project/DEG_analysis_disease/results/Dis.vs.all/DE_results/CELLTYPE/'

chris.genes <- fread('/home/USSR/awc30/liver_project/manuscript_figures/DEGs_for_figure_S2.csv')

cellTypes <- list('B-cell1'='B-cell.1', 
                  'B-cell2'='B-cell.2', 
                  'Chol'='Cholangiocytes', 
                  'Endo'='Endothelial', 
                  'Lymph'='Lymphocytes', 
                  'Macro'='Macrophages', 
                  'Neutro'='Neutrophils',
                  'Stell'='Stellate')

celltype <- names(cellTypes)[3]
for (celltype in names(cellTypes)) {
  print(celltype)
  
  seuFile <- paste0('Aggr_Jan2022_', celltype, '_harmony_th=0.rds')

  seu <- readRDS(paste0(projDir, seuDir, seuFile))

  disease.order <- c('Healthy control', 'NAFLD', 'NASH w/o cirrhosis',
                     'NASH with cirrhosis', 'end stage')
  seu@meta.data$Disease.status <- factor(seu@meta.data$Disease.status,
                                         levels=disease.order)
  m <- seu@meta.data
  
  # supp 2A type
  L <- make_UMAP_plot(seu, 
                 'umap_harmony_t.0', 
                 'Disease.status')
  
  L[[1]]
  pdf(file=paste0(plotdir, 'S2A_', celltype, '.pdf'))
    print(L[[2]])
    print(L[[1]])
  dev.off()
  
  
  # 2B : DEG bubble plot and 3 DEG in each disease status
  DE.dir <- '/home/USSR/awc30/liver_project/DEG_analysis_disease/results/Dis.vs.all/DE_results/'
  cellDir <- cellTypes[[celltype]]
  dis.comparisons <- c('Healthy.control',
                       'NAFLD',
                       'NASH.w.o.cirrhosis',
                       'NASH.with.cirrhosis',
                       'end.stage')
  
  dis.state <- dis.comparisons[4]
  plt.genes <- c()
  # for (dis.state in dis.comparisons) {
  #   DE.table.path <- paste0(DE.dir, cellDir, '/all_vs_', dis.state, '_DE_wilcox.csv')
  #   DE.table <- fread(DE.table.path)
  #   DE.table <- DE.table[order(DE.table$avg_log2FC, decreasing=T)]
  #   if (nrow(DE.table) > 0) {
  #     curr.up.genes <- DE.table$gene[1:min(3, nrow(DE.table))]
  #     plt.genes <- c(plt.genes, curr.up.genes)
  #   }
  # }
  # plot the genes Chris requested
  plt.genes <- chris.genes$gene[chris.genes$cellName==celltype]
  
  plt.genes <- unique(plt.genes)
  Idents(seu) <- seu@meta.data$Disease.status
  p <- DotPlot(seu, assay='RNA', features=plt.genes,
               scale=T)
  p <- p+coord_flip()+
    xlab('')+
    ylab('')+
    theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
  p
  p.nolab <- p + theme(legend.position = 'none')

  pdf(file=paste0(plotdir, 'S2B_', celltype, '.pdf'), height=7, width=5)
    print(p)
    print(p.nolab)
  dev.off()
}









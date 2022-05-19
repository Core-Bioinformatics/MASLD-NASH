rm(list=ls())
library(cowplot)
library(Seurat)
library(ggplot2)
library(data.table)
source('/home/USSR/awc30/liver_project/manuscript_figures/manuscript_UMAPs&dotplots/plot_functions.R')


# Figure 1 & supp. figure 1 - umaps and dotplots using all cells ---------------

seu.data.dir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'
seu.data.version.dir <- '6.Aggr-Jan2022/'
seu.filename <- 'Aggr_Jan2022_AllCells_harmony_th=0.rds'
plotdir <- '/home/USSR/awc30/liver_project/manuscript_figures/manuscript_UMAPs&dotplots/plots/'

seu <- readRDS(paste0(seu.data.dir, seu.data.version.dir, seu.filename))
dim(seu@assays$RNA@counts)

length(unique(seu@meta.data$orig.ident))

m <- seu@meta.data
harmony_reduction <- 'umap_harmony_t.0'
table(seu@meta.data$cell.annotation)

# remove the 122 unknown celltype cells
Idents(seu) <- m$cell.annotation
seu.filt <- subset(seu, idents='unknown', invert=T)


# sort out some factor levels
seu.filt@meta.data$Fibrosis.score..F0.4. <- as.factor(seu.filt@meta.data$Fibrosis.score..F0.4.)
disease.order <- c('Healthy control', 'NAFLD', 'NASH w/o cirrhosis',
                   'NASH with cirrhosis', 'end stage')
seu.filt@meta.data$Disease.status <- factor(seu.filt@meta.data$Disease.status,
                                            levels=disease.order)





# The number of nuclei recovered per sample ------
m <- seu.filt@meta.data
m$barcode <- rownames(m)
m <- data.table(m)

nCells <- m[, .('nCells'=.N), by=.(orig.ident, Disease.status)]

# fix sample names to get rid of SLX-, so consistent with GEO
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19591-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19693-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19750-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19940-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20150-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20266-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20270-', 'B-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20289-', 'C-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20290-', 'D-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20793-', 'E-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20985-', 'F-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-21151-', 'G-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-21153-', 'H-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-21155-', 'I-')


nCells$Disease.status <- factor(nCells$Disease.status, 
                                levels=c('Healthy control', 'NAFLD', 'NASH with cirrhosis', 
                                'NASH w/o cirrhosis', 'end stage'))
p <- ggplot(nCells, aes(x=Disease.status, y=nCells, color=Disease.status))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(size=1, alpha=0.6)+
  theme_bw()+
  ylab('num. nuclei')+
  xlab('')+
  ggtitle(paste0('Overall mean is ', round(mean(nCells$nCells)), ' nuclei'))+
  theme(legend.position='none', 
        axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p
ggsave(filename=paste0(plotdir, 'num_nuclei_by_sample.pdf'),
       width=4, height=4)


# the number of genes expressed per cell per sample --------------
m <- seu.filt@meta.data
m$barcode <- rownames(m)
m <- data.table(m)

nCells <- m[, .('nGenes.per.cell'=mean(nFeature_RNA)), by=.(orig.ident, Disease.status)]

# fix sample names to get rid of SLX-, so consistent with GEO
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19591-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19693-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19750-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-19940-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20150-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20266-', '')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20270-', 'B-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20289-', 'C-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20290-', 'D-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20793-', 'E-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-20985-', 'F-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-21151-', 'G-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-21153-', 'H-')
nCells$orig.ident <- str_replace_all(nCells$orig.ident, 'SLX-21155-', 'I-')


nCells$Disease.status <- factor(nCells$Disease.status, 
                                levels=c('Healthy control', 'NAFLD', 'NASH with cirrhosis', 
                                         'NASH w/o cirrhosis', 'end stage'))
p <- ggplot(nCells, aes(x=Disease.status, y=nGenes.per.cell, color=Disease.status))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(size=1, alpha=0.6)+
  theme_bw()+
  ylab('mean genes expressed per nuclei')+
  xlab('')+
  ggtitle(paste0('Overall mean is: ', round(mean(nCells$nGenes.per.cell)), ' genes'))+
  theme(legend.position='none', 
        axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p
ggsave(filename=paste0(plotdir, 'num_genes_per_nuclei_by_sample.pdf'),
       width=4, height=4)



# 1D -----

# manually set colours for the cells to make the small groups pop
cols <- gg_color_hue(9)
cellType.cols <- c(cols[9], cols[5], cols[3], cols[1], cols[6], cols[7],
              cols[8], cols[4], cols[2])
L <- make_UMAP_plot(seu.filt, 
                    'umap_harmony_t.0', 
                    'cell.annotation', 
                    cellType.cols)
L[[1]]
pdf(file=paste0(plotdir, 'fig1D.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()


# 1F ------
L <- make_UMAP_plot(seu.filt, 
                    'umap_harmony_t.0', 
                    'Patient.ID')
L[[1]]
pdf(file=paste0(plotdir, 'fig1F.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()

# alternatively, facet by Disease state as well
umap.reduction='umap_harmony_t.0'
group.by.col='Patient.ID'
facet.col='Disease.status'

L <- make_facet_umap(seu, umap.reduction, 
                     group.by.col, facet.col, 
                     nrow=1)
L[[1]]
pdf(file=paste0(plotdir, 'fig1F_alternative.pdf'), 
    width=25, height=5)
    print(L[[2]])
    print(L[[1]])
dev.off()




# 1G ------
L <- make_UMAP_plot(seu.filt, 
                    'umap_harmony_t.0', 
                    'Disease.status')
L[[1]]
pdf(file=paste0(plotdir, 'fig1G.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()


# 1E -------
marker.genes <- c('IGKC', 'MZB1', 'BANK1', 'CD247', 'CD2', 'COL3A1', 'DCN','MNDA', 
                  'FCN1', 'CD163', 'MARCO', 'CFTR', 'KRT7', 'STAB2', 'ASGR1', 'CYP3A4')

# celltype.order <- c('Hepatocytes', 'Endothelial','Cholangiocytes', 
#                     'Macrophagees','Neutrophils', 'Stellate', 'Lymphocytes', 
#                     'B-cell 1', 'B-cell 2')

# make it without the unknown cells
p <- DotPlot(seu.filt, assay='SCT',
             features=marker.genes, 
             scale=T)
p <- p+coord_flip()+
  xlab('')+
  ylab('')+
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p.nolab <- p + theme(legend.position = 'none')

pdf(file=paste0(plotdir, 'fig1E.pdf'), height=5)
  print(p)
  print(p.nolab)
dev.off()

# make it with the 112 unknown cells
p <- DotPlot(seu, assay='SCT',
             features=marker.genes, 
             scale=T)
p <- p+coord_flip()+
  xlab('')+
  ylab('')+
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
p.nolab <- p + theme(legend.position = 'none')

pdf(file=paste0(plotdir, 'fig1E_w_unknown.pdf'), height=5)
  print(p)
  print(p.nolab)
dev.off()


# Supplementary for fig 1 ------------------------------------------------------

# Uncorrected UMAPs ------------------

# cell annotation --------
L <- make_UMAP_plot(seu.filt, 
                    'umap', 
                    'cell.annotation',
                    cellType.cols)
L[[2]]
pdf(file=paste0(plotdir, 'S1_cellAnnotation.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()

# patient --------
L <- make_UMAP_plot(seu.filt, 
                    'umap', 
                    'Patient.ID')
L[[2]]
pdf(file=paste0(plotdir, 'S1_patient.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()

# alternatively, facet by Disease state as well
umap.reduction='umap_harmony_t.0'
group.by.col='Patient.ID'
facet.col='Disease.status'

L <- make_facet_umap(seu, 'umap', 
                     group.by.col, facet.col, 
                     nrow=1)
L[[1]]
pdf(file=paste0(plotdir, 'S1_patient_alternative.pdf'), 
    width=25, height=5)
print(L[[2]])
print(L[[1]])
dev.off()

# disease state ---------
L <- make_UMAP_plot(seu.filt, 
                    'umap', 
                    'Disease.status')
L[[2]]
pdf(file=paste0(plotdir, 'S1_diseaseStatus.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()

# umap by fibrosis score ------
L <- make_UMAP_plot(seu.filt, 
                    'umap_harmony_t.0',
                    'Fibrosis.score..F0.4.')
L[[1]] <- L[[1]]+scale_color_viridis_d()
L[[2]] <- L[[2]]+scale_color_viridis_d()

pdf(file=paste0(plotdir, 'S1_fibrosis_score.pdf'))
  print(L[[2]])
  print(L[[1]])
dev.off()



# facetted disease stats -------
p <- DimPlot(seu.filt, 
             reduction='umap_harmony_t.0',
             group.by = 'Disease.status',
             split.by='Disease.status', 
             ncol=3)
p.noleg <- p+theme(legend.position='none')

pdf(file=paste0(plotdir, 'S1_DisStatus_facetted.pdf'), height=5.5, width=7)
  print(p.noleg)
  print(p)
dev.off()


# UMAPs with expression of the markers from 1E
marker.genes <- c('IGKC', 'MZB1', 'BANK1', 'CD247', 'CD2', 'COL3A1', 'DCN','MNDA', 
                  'FCN1', 'CD163', 'MARCO', 'CFTR', 'KRT7', 'STAB2', 'ASGR1', 'CYP3A4')

plist <- list()
mg <- marker.genes[3]
for (mg in marker.genes) {
  p <- FeaturePlot(seu.filt, reduction='umap_harmony_t.0', features=mg, slot='data')
  # p <- p+scale_color_viridis_c()+
  p <- p+ xlab('UMAP1') +
    ylab('UMAP2')+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank())
  plist[[mg]] <- p
}

pdf(file=paste0(plotdir, 'S1_Marker_UMAPs.pdf'), height=7, width=7)
for (mg in names(plist)) {
    print(plist[[mg]])
  }
dev.off()






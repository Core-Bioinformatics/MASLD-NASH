rm(list=ls())

library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)
library(lme4)
library(cowplot)

source('/home/USSR/awc30/liver_project/manuscript_figures/manuscript_UMAPs&dotplots/plot_functions.R')
source('/home/USSR/awc30/liver_project/DEG_analysis_clusters/scripts/do_DGE_functions.R')



projDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'
seuDir <- '/6.Aggr-Jan2022/'
plotDir <- '/home/USSR/awc30/liver_project/manuscript_figures/manuscript_UMAPs&dotplots/plots/fig3/'

dir.create(plotDir)


# see powerpoint figure 3 explanations for what Chris has requested for this
# "slides" refer to this powerpoint ---


# Hypoxis/oxstress marker genes -----------------------------------------------
hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))

gene.vec <- c('HIF1A','ARNT','IGFBP1','BNIP3','PDK1','SLC2A1','ALDOA','ENO1', 'LDHA',
              'AOX1','GST1A','SOD1','SOD2','GPX3','TXN','TXNRD1','GAD1','DUSP1',
              'HSPA4','HSPA2','HSPA5','HSPA8','HSPA9','HSF1',
              'VIM','SNAI1','SNAI2','CDH1','CDH2','TWIST1','ZEB2','ZEB1',
              'CDKN2A','CDKN2B','CDKN1A','CDKN1B','TP53',
              'MKI67','CDK1','PCNA','CCNE1','MCM3',
              'CXCL8','IL6','TNF','TNFRSF1A','TNFRSF1B','CRP','CXCL1','CCL20','LIF',
              'YAP','TAZ','BCL2','CCN1','CCN2','E2F1','CCND1','AREG','MCL1','FGF1')
gene.vec <- sort(unique(gene.vec))

gv <- gene.vec[1]
pdf(paste0(plotDir, 'stress_hypoxia_oxstress_UMAP_hep.pdf'), w=10, h=10)
for (gv in gene.vec) {
  if (gv %in% rownames(hep.seu@assays$SCT@data)) {
    p <- make_gene_expression_umap(hep.seu, gv, is.log=T)
    print(p)
  }
}
dev.off()

pdf(paste0(plotDir, 'stress_hypoxia_oxstress_UMAP_chol.pdf'), w=10, h=10)
for (gv in gene.vec) {
  if (gv %in% rownames(chol.seu@assays$SCT@data)) {
    p <- make_gene_expression_umap(chol.seu, gv, is.log=T)
    print(p)
  } 
}
dev.off()



# Stem / progenitor marker genes -----------------------------------------------
hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))

gene.vec <- c('MKI67', 'AFP', 'DLK1', 'TBX3', 'NANOG', 'POU5F1', 'PROM1', 'SPINK1', 
              'CAV1', 'GPRC5B', 'MCAM', 'NCAM1', 'STAT4', 
              'TACSTD2', 'EPCAM', 'PDX1', 'LGR5', 'PROM1', 'SPINK1', 'SPP1', 'PROX1',
              'AFP', 'DLK1', 'FOXM1', 'DONSON', 'NPW', 'SALL4', 'COL4A2', 'CBX1', 'FST', 
              'TERT')
gene.vec <- sort(unique(gene.vec))

gv <- gene.vec[1]
pdf(paste0(plotDir, 'stem_progenitor_marker_UMAP_hep.pdf'), w=10, h=10)
for (gv in gene.vec) {
  if (gv %in% rownames(hep.seu@assays$SCT@data)) {
    p <- make_gene_expression_umap(hep.seu, gv, is.log=T)
    print(p)
  }
}
dev.off()

pdf(paste0(plotDir, 'stem_progenitor_marker_UMAP_chol.pdf'), w=10, h=10)
for (gv in gene.vec) {
  if (gv %in% rownames(chol.seu@assays$SCT@data)) {
    p <- make_gene_expression_umap(chol.seu, gv, is.log=T)
    print(p)
  } 
}
dev.off()


# oxidative stress, hypoxia, stress and epithelial-to-mesenchymal transition (EMT) ------
hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))

gene.vec <- c('HIF1A', 'ARNT', 'IGFBP1', 'BNIP3', 'PDK1', 'SLC2A1', 'ALDOA', 'ENO1', 'LDHA',
              'AOX1', 'GST1A', 'SOD1', 'SOD2', 'GPX3', 'TXN', 'TXNRD1', 'GAD1', 'DUSP1',
              'HSPA4', 'HSPA2', 'HSPA5', 'HSPA8', 'HSPA9', 'HSF1', 
              'VIM', 'SNAI1', 'SNAI2', 'CDH1', 'CDH2', 'TWIST1', 'ZEB2', 'ZEB1'
              )
gene.vec <- sort(unique(gene.vec))

gv <- gene.vec[1]
pdf(paste0(plotDir, 'stress_marker_UMAP_hep.pdf'), w=10, h=10)
for (gv in gene.vec) {
  if (gv %in% rownames(hep.seu@assays$SCT@data)) {
    p <- make_gene_expression_umap(hep.seu, gv, is.log=T)
    print(p)
  }
}
dev.off()

pdf(paste0(plotDir, 'stress_marker_UMAP_chol.pdf'), w=10, h=10)
for (gv in gene.vec) {
  if (gv %in% rownames(chol.seu@assays$SCT@data)) {
    p <- make_gene_expression_umap(chol.seu, gv, is.log=T)
    print(p)
  } 
}
dev.off()

# Slide 2 ----------------------------------------------------------------------
hep.seu.file <- 'Aggr_Jan2022_Hep_harmony_th=0.rds'
hep.seu <- readRDS(paste0(projDir, seuDir, hep.seu.file))
m <- hep.seu@meta.data

# UMAP
p <- DimPlot(hep.seu, reduction='umap_harmony_t.0', group.by='SCT_snn_harmony_t.0.0.4')
p
ggsave(paste0(plotDir, 'slide_2-1.pdf'), width=10, height=10)

# tmp <- FeaturePlot(hep.seu, reduction='umap_harmony_t.0', features=c('KRT7'))
# tmp

# Heatmap
gene.vec <- c("KRT7", "CFTR", "EPCAM", "AQP1", "ASGR1", "TTR", "CYP3A4", "HNF4A" )
p <- make_expression_heatmap(hep.seu, gene.vec, 'SCT_snn_harmony_t.0.0.4',
                             scale.by.gene=T)
p
ggsave(paste0(plotDir, 'slide_2-2.pdf'))

# Heatmap with alternative genes
gene.vec <- c("ALB", "ASGR1", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "CFTR", "EPCAM", "AQP1")
p <- make_expression_heatmap(hep.seu, gene.vec, 'SCT_snn_harmony_t.0.0.4',
                             scale.by.gene=T)
p
ggsave(paste0(plotDir, 'slide_2-2_alternative.pdf'))

rm(hep.seu, m)


# expression UMAP for each gene in the heatmap
gene.vec <- c("ALB", "ASGR1", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "CFTR", "EPCAM", "AQP1")

gv <- gene.vec[1]
pdf(paste0(plotDir, 'slide_2-3.pdf'), w=10, h=10)
for (gv in gene.vec) {
  p <- make_gene_expression_umap(hep.seu, gv, is.log=T)
  print(p)
}
dev.off()


# Slide 3 ----------------------------------------------------------------------
seu.file <- 'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'
seu <- readRDS(paste0(projDir, seuDir, seu.file))
m <- seu@meta.data

p <- DimPlot(seu, reduction='umap_harmony_t.0', group.by='clust.9.subcluster',
             pt.size=0.8, order=c('0', '1', '2', '3')) +
  scale_color_manual(values=c('0'='dodgerblue4',
                              '1'='dodgerblue3',
                              '2'='deepskyblue2',
                              '3'='deepskyblue',
                              'NO'='grey87'))+
  theme(legend.position='bottom',
        legend.text=element_text(size=6))
p
ggsave(paste0(plotDir, 'slide_3-1.pdf'), width=5, height=6)


gene.vec <- c("KRT7", "CFTR", "EPCAM", "AQP1", "ASGR1", "TTR", "CYP3A4", "HNF4A" )
p <- make_expression_heatmap(seu, gene.vec, 'clust.9.subcluster',
                             scale.by.gene=T)
p
ggsave(paste0(plotDir, 'slide_3-2.pdf'))

# Heatmap with alternative genes
gene.vec <- c("ALB", "ASGR1", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "CFTR", "EPCAM", "AQP1")
p <- make_expression_heatmap(seu, gene.vec, 'clust.9.subcluster',
                             scale.by.gene=T)
p
ggsave(paste0(plotDir, 'slide_3-2_alternative.pdf'))


p <- make_gene_expression_umap(seu, 'KRT7')
p
ggsave(paste0(plotDir, 'slide_3-3.pdf'), width=10, height=10)

# Slide 4 ----------------------------------------------------------------------
seu.file <- 'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'
seu <- readRDS(paste0(projDir, seuDir, seu.file))

# Hepatocytes in cluster 9 subclusters ----
hep.m <- data.table(seu@meta.data[, c('Patient.ID', 'Disease.status', 
                                          'clust.9.subcluster')])
hep.m$new.name <- 'hepatocytes'
hep.m$new.name[hep.m$clust.9.subcluster %in% c('2', '3')] <- 'chol. like hep.'

clust.id <- 'chol. like hep.'

cluster.prop.dt <- get_proportion_of_cells_in_cluster(hep.m, 
                                                      'new.name',
                                                      clust.id)
p.bar <- make_barplot(cluster.prop.dt, '')
p.bar <- p.bar +
  ylim(0, 0.03) # add a bit of room so *s can be added manually
p.bar

pVal.df <- calc_all_pairwise_disease_pvals(hep.m, 'new.name', clust.id)
table(hep.m[, c('Patient.ID', 'new.name')])
p.prob.plot <- plot_probability(pVal.df)
p.prob.plot

p.both <- plot_grid(p.bar, p.prob.plot, nrow=1)

ggsave(filename=paste0(plotDir, 'slide_4.pdf'),
       width=10, height=4)


# Hepatocytes in cluster c16 (@1.6 res) subcluster ----
hep.m <- data.table(seu@meta.data[, c('Patient.ID', 'Disease.status', 
                                      'SCT_snn_harmony_t.0.1.6')])

hep.m$new.name <- 'hepatocytes'
hep.m$new.name[hep.m$SCT_snn_harmony_t.0.1.6 %in% c('16')] <- 'c16'
unique(hep.m$new.name)
table(hep.m[, c('Disease.status', 'new.name')])


clust.id <- 'c16'
col.id <- 'new.name'
cluster.prop.dt <- get_proportion_of_cells_in_cluster(hep.m, 
                                                      col.id,
                                                      clust.id)
p.bar <- make_barplot(cluster.prop.dt, '')
p.bar
p.bar <- p.bar +
  ylim(0, 0.075) # add a bit of room so *s can be added manually - CAREFUL - this can remove bars if too low!
p.bar

pVal.df <- calc_all_pairwise_disease_pvals(hep.m, 'new.name', clust.id)
table(hep.m[, c('Patient.ID', 'new.name')])
p.prob.plot <- plot_probability(pVal.df)
p.prob.plot

p.both <- plot_grid(p.bar, p.prob.plot, nrow=1)

ggsave(filename=paste0(plotDir, 'slide_4_c16.pdf'),
       width=10, height=4)


# Hepatocytes in c16 + c9.2 + c9.3 subcluster ----

hep.m <- data.table(seu@meta.data[, c('Patient.ID', 'Disease.status', 
                                      'SCT_snn_harmony_t.0.1.6', "clust.9.subcluster")])

clust.id <- 'chol. like hep.'
hep.m$new.name <- 'hepatocytes'
hep.m$new.name[hep.m$SCT_snn_harmony_t.0.1.6 %in% c('16')] <- clust.id
hep.m$new.name[hep.m$clust.9.subcluster %in% c('2', '3')] <- clust.id
unique(hep.m$new.name)
table(hep.m[, c('Disease.status', 'new.name')])


col.id <- 'new.name'
cluster.prop.dt <- get_proportion_of_cells_in_cluster(hep.m, 
                                                      col.id,
                                                      clust.id)
p.bar <- make_barplot(cluster.prop.dt, '')
p.bar
p.bar <- p.bar +
  ylim(0, 0.1) # add a bit of room so *s can be added manually - CAREFUL - this can remove bars if too low!
p.bar

pVal.df <- calc_all_pairwise_disease_pvals(hep.m, 'new.name', clust.id)
table(hep.m[, c('Patient.ID', 'new.name')])
p.prob.plot <- plot_probability(pVal.df)
p.prob.plot

p.both <- plot_grid(p.bar, p.prob.plot, nrow=1)

ggsave(filename=paste0(plotDir, 'slide_4_c9.2,c9.3,c16.pdf'),
       width=10, height=4)


# Slide 5 ----------------------------------------------------------------------

# do DGE of hep clusters 9.2, 9.3 combined vs all heps -------------------------
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'))

# make combined 9.2, 9.3, and all others are Hep
seu@meta.data$hep.cluster <- 'hepatocytes'
seu@meta.data$hep.cluster[seu@meta.data$clust.9.subcluster %in% c('2', '3')] <- 'chol. like hep.'
unique(seu@meta.data$hep.cluster)

Idents(seu) <- seu@meta.data$hep.cluster
# DE parameters:
logFC.th=0.5
p.val.adj.th=0.1

all.markers <- FindAllMarkers(seu,
                              logfc.threshold=logFC.th,
                              test.use='wilcox',
                              min.pct=0.1,
                              return.thres=p.val.adj.th,
                              base=2) # base = which to use for logs
fwrite(file=paste0(plotDir, 'slide_5_DE_table.csv'), 
       all.markers)
all.markers <- fread(file=paste0(plotDir, 'slide_5_DE_table.csv'))

gene.vec <- c("BICC1", "DCDC2", "FGF13", "FMNL2", "FGFR2", "SERPINE1", 'KLF6', "BCL2", "CDH6", "CREB5", "KRT23", 'SOX4')
tmp <- all.markers[all.markers$gene %in% gene.vec,]

p <- make_expression_heatmap(seu, gene.vec, 'hep.cluster',
                             scale.by.gene=F, 
                             ) # don't scale with only two columns!
# compress colorbar so can see other gene differences for all genes
p <- p+scale_fill_viridis_c(limits=c(0, 4), oob=scales::squish,
                            breaks=c(0,1,2,3,4),
                            labels=c('0','1','2','3','> 4'))
p
ggsave(paste0(plotDir, 'slide_5-1+SOX4+KLF6.pdf'), width=4)


# expression UMAP for each gene in the heatmap
seu.hep <- seu
gene.vec <- c("BICC1", "DCDC2", "FGF13", "FMNL2", "FGFR2", "SERPINE1", 'KLF6', "BCL2", "CDH6", "CREB5", "KRT23", 'SOX4')
gv <- gene.vec[1]
pdf(paste0(plotDir, 'slide_5-2+SOX4+KLF6.pdf'), w=10, h=10)
for (gv in gene.vec) {
  p <- make_gene_expression_umap(seu.hep, gv, is.log=T)
  print(p)
}
dev.off()

seu.chol <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))
gene.vec <- c("BICC1", "DCDC2", "FGF13", "FMNL2", "FGFR2", "SERPINE1", 'KLF6', "BCL2", "CDH6", "CREB5", "KRT23", 'SOX4')
gv <- gene.vec[1]
pdf(paste0(plotDir, 'slide_5-3+SOX4+KLF6.pdf'), w=10, h=10)
for (gv in gene.vec) {
  p <- make_gene_expression_umap(seu.chol, gv, is.log=T)
  print(p)
}
dev.off()


# Slide 6 ----------------------------------------------------------------------
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m <- seu@meta.data
p <- DimPlot(seu, reduction='umap_harmony_t.0', group.by='SCT_snn_harmony_t.0.0.4')
p
ggsave(filename=paste0(plotDir, 'slide_6-1.pdf'), width=10, height=10)

# heatmap of expression
gene.vec <- c("ASGR1", "ALB", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "SOX9", "CFTR")
p <- make_expression_heatmap(seu, gene.vec, 'SCT_snn_harmony_t.0.0.4',
                             scale.by.gene=T) 
p
ggsave(paste0(plotDir, 'slide_6-2.pdf'))

p <- p+scale_fill_viridis_c(limits=c(-1, 1), oob=scales::squish,
                            breaks=c(-1, 0, 1),
                            labels=c('< -1', '0', '> 1'))
p
ggsave(paste0(plotDir, 'slide_6-2_mod_scale.pdf'))

# same heatmap with alternative genes
gene.vec <- c("ALB", "ASGR1", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "CFTR", "EPCAM", "AQP1")
p <- make_expression_heatmap(seu, gene.vec, 'SCT_snn_harmony_t.0.0.4',
                             scale.by.gene=T) 
p
ggsave(paste0(plotDir, 'slide_6-2_alternative.pdf'))
p <- p+scale_fill_viridis_c(limits=c(-1, 1), oob=scales::squish,
                            breaks=c(-1, 0, 1),
                            labels=c('< -1', '0', '> 1'))
p
ggsave(paste0(plotDir, 'slide_6-2_mod_scale_alternative.pdf'))


p <- make_gene_expression_umap(seu, 'ASGR1', is.log=F)
p
ggsave(paste0(plotDir, 'slide_6-3.pdf'), width=10, height=10)


p <- make_gene_expression_umap(seu, 'GPC6', is.log=F)
p
ggsave(paste0(plotDir, 'slide_6-4.pdf'), width=10, height=10)

# expression UMAP for each gene in the heatmap
gene.vec <- c("ALB", "ASGR1", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "CFTR", "EPCAM", "AQP1")

gv <- gene.vec[1]
pdf(paste0(plotDir, 'slide_6-5.pdf'), w=10, h=10)
for (gv in gene.vec) {
  p <- make_gene_expression_umap(seu, gv, is.log=T)
  print(p)
}
dev.off()


# Slide 7 ----------------------------------------------------------------------
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m <- seu@meta.data

# interested in clusters at 0.4 res and 1.6 res ----
chol.m <- data.table(seu@meta.data[, c('Patient.ID', 'Disease.status', 
                                      'SCT_snn_harmony_t.0.0.4', 
                                      'SCT_snn_harmony_t.0.1.6')])

table(chol.m[, c('Disease.status', 'SCT_snn_harmony_t.0.0.4')])


clust.ids <- list('1'='SCT_snn_harmony_t.0.0.4', '5'='SCT_snn_harmony_t.0.0.4',
                  '9'='SCT_snn_harmony_t.0.0.4',
                  '22'='SCT_snn_harmony_t.0.1.6')

plt.list <- list()
for (clust.id in names(clust.ids)) {
  curr.col <- clust.ids[[clust.id]]
  
  cluster.prop.dt <- get_proportion_of_cells_in_cluster(chol.m, 
                                                        curr.col,
                                                        clust.id)
  
  p.bar <- make_barplot(cluster.prop.dt, paste0('Chol. cells in cluster ', clust.id))
  p.bar
  p.bar <- p.bar +
    scale_y_continuous(expand=expansion(mult=c(.05,.4))) # add a bit of room so *s can be added manually
  p.bar
  pVal.df <- calc_all_pairwise_disease_pvals(chol.m, curr.col, clust.id)
  #table(hep.m[, c('Patient.ID', 'new.name')])
  p.prob.plot <- plot_probability(pVal.df)
  p.prob.plot
  
  p.both <- plot_grid(p.bar, p.prob.plot, nrow=1)
  plt.list[[clust.id]] <- p.both 
}

pdf(file=paste0(plotDir, 'slide_7.pdf'), w=10, h=4)
  for (n in names(plt.list)) {
    p <- plt.list[[n]]
    print(p)
  } 
dev.off()


# Cholangiocytes in combined c5, c9 (0.4 res), c22, c3, c7 (1.6 res) subclusters ----
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m <- seu@meta.data

chol.m <- data.table(seu@meta.data[, c('Patient.ID', 'Disease.status', 
                                      'SCT_snn_harmony_t.0.1.6', "SCT_snn_harmony_t.0.0.4")])

clust.id <- 'hep. like chol.'
chol.m$new.name <- 'cholangiocytes'
chol.m$new.name[chol.m$SCT_snn_harmony_t.0.1.6 %in% c('22', '3', '7')] <- clust.id
chol.m$new.name[chol.m$SCT_snn_harmony_t.0.0.4 %in% c('5', '9')] <- clust.id
unique(chol.m$new.name)
table(chol.m[, c('Disease.status', 'new.name')])
table(chol.m[, c('SCT_snn_harmony_t.0.1.6', 'new.name')])
table(chol.m[, c('SCT_snn_harmony_t.0.0.4', 'new.name')])


col.id <- 'new.name'
cluster.prop.dt <- get_proportion_of_cells_in_cluster(chol.m, 
                                                      col.id,
                                                      clust.id)
p.bar <- make_barplot(cluster.prop.dt, '')
p.bar
p.bar <- p.bar +
  ylim(0, 0.35) # add a bit of room so *s can be added manually - CAREFUL - this can remove bars if too low!
p.bar

pVal.df <- calc_all_pairwise_disease_pvals(chol.m, 'new.name', clust.id)
table(chol.m[, c('Patient.ID', 'new.name')])
p.prob.plot <- plot_probability(pVal.df)
p.prob.plot

p.both <- plot_grid(p.bar, p.prob.plot, nrow=1)
p.both

ggsave(filename=paste0(plotDir, 'slide_7_c5,c9,c22,c3,c7.pdf'),
       width=10, height=4)




# Slide 8 ----------------------------------------------------------------------

# Boxplot of Chol and hep marker gene expression in cholangiocytes, hepatocytes,
# and subclusters --------------------------------------------------------------
hep.seu <- readRDS(file=paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0_c5_subclusters.rds'))
chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'))

# add subcluster info ---
chol.hep.seu@meta.data$subcluster <- chol.hep.seu@meta.data$cell.annotation

# add hep subcluster info
subclust='0'
for (subclust in c('2', '3')) {
    barcodes <- rownames(hep.seu@meta.data[hep.seu@meta.data$clust.9.subcluster==subclust,])
    chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- 'chol. like heps'
}
unique(chol.hep.seu@meta.data$subcluster)

# add chol 5 subcluster info
for (subclust in unique(chol.seu@meta.data$cluster.5.subcluster)) {
  if (subclust=='NO') {
    next
  } else {
    barcodes <- rownames(chol.seu@meta.data[chol.seu@meta.data$cluster.5.subcluster==subclust,])
    chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- paste0('chol-5')
  }
}

# plot with and without hepatocytes ---
ordered.markers <- c('KRT7', 'CFTR', 'CTNND2', 'ALB', 'ASGR1', 'TTR')
ordered.subclusters <- c('Cholangiocytes', 'chol-5',
                         'chol. like heps', 'Hepatocytes')
p <- plot_violinPlot(chol.hep.seu, ordered.markers, ordered.subclusters)
p
ggsave(paste0(plotDir, 'slide_8_w_heps.pdf'), 
       width=7, height=5)

chol.hep.seu.sub <- subset(chol.hep.seu, subset=subcluster=='Hepatocytes', invert=T)
p <- plot_violinPlot(chol.hep.seu.sub, ordered.markers, c('Cholangiocytes', 'chol-5', 'chol. like heps'))
ggsave(paste0(plotDir, 'slide_8_wo_heps.pdf'), 
       width=7, height=5)


# add chol subcluster 22 info as well ---
barcodes <- rownames(chol.seu@meta.data[chol.seu@meta.data$SCT_snn_harmony_t.0.1.6=="22",])
chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- paste0('chol-22')

ordered.markers <- c('KRT7', 'CFTR', 'CTNND2', 'ALB', 'ASGR1', 'TTR')
ordered.subclusters <- c('Cholangiocytes', 'chol-5', 'chol-22',
                         'chol. like heps', 'Hepatocytes')
p <- plot_violinPlot(chol.hep.seu, ordered.markers, ordered.subclusters)
p
ggsave(paste0(plotDir, 'slide_8_w_heps+chol22.pdf'), 
       width=7, height=5)

chol.hep.seu.sub <- subset(chol.hep.seu, subset=subcluster=='Hepatocytes', invert=T)
p <- plot_violinPlot(chol.hep.seu.sub, ordered.markers, c('Cholangiocytes', 'chol-5', 'chol-22', 'chol. like heps'))
ggsave(paste0(plotDir, 'slide_8_wo_heps+chol22.pdf'), 
       width=7, height=5)



# slide 9 ----------------------------------------------------------------------
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m <- seu@meta.data
p <- DimPlot(seu, reduction='umap_harmony_t.0', group.by='SCT_snn_harmony_t.0.1.6')
p
ggsave(filename=paste0(plotDir, 'slide_9-1.pdf'), width=10, height=10)

gene.vec <- c("ASGR1", "ALB", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "SOX9", "CFTR")
seu.sub <- subset(seu, subset=SCT_snn_harmony_t.0.0.4=='1')

p <- make_expression_heatmap(seu.sub, gene.vec, 'SCT_snn_harmony_t.0.1.6',
                             scale.by.gene=T) 
p
ggsave(paste0(plotDir, 'slide_9-2.pdf'))

gene.vec <- c("ALB", "ASGR1", "TTR", "ASS1", "PCK1", "ABCC2", "GPC6", "HNF4A", 
              "KRT19", "KRT7", "CFTR", "EPCAM", "AQP1")

p <- make_expression_heatmap(seu.sub, gene.vec, 'SCT_snn_harmony_t.0.1.6',
                             scale.by.gene=T) 
p
ggsave(paste0(plotDir, 'slide_9-2_alternative.pdf'))




# slide 10 ---------------------------------------------------------------------
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))
seu@meta.data$is.clust.22 <- seu@meta.data$cell.annotation
seu@meta.data$is.clust.22[seu@meta.data$SCT_snn_harmony_t.0.1.6=='22'] <- 'clust-22'
table(seu@meta.data[, c('is.clust.22', 'SCT_snn_harmony_t.0.1.6')])

Idents(seu) <- seu@meta.data$is.clust.22
# DE parameters:
logFC.th=0.5
p.val.adj.th=0.1

all.markers <- FindMarkers(seu,
                           ident.1='clust-22',
                           ident.2='Cholangiocytes',
                           logfc.threshold=logFC.th,
                           test.use='wilcox',
                           min.pct=0.1,
                           return.thres=p.val.adj.th,
                           base=2) # base = which to use for logs
all.markers$comparison='clust-22 vs cholangiocytes'
all.markers$gene <- rownames(all.markers)
all.markers[['+ve logFC means greater in']] <- 'clust-22'
fwrite(file=paste0(plotDir, 'slide_10_DE_table.csv'),
       all.markers)
all.markers <- fread(file=paste0(plotDir, 'slide_10_DE_table.csv'))

# gene.vec <- c("BICC1", "DCDC2", "FGF13", "FMNL2", "FGFR2", "SERPINE1", "BCL2", "CDH6", "CREB5", "KRT23", 'SOX4')
# tmp <- all.markers[all.markers$gene %in% gene.vec,]
# 
# p <- make_expression_heatmap(seu, gene.vec, 'is.clust.22',
#                              scale.by.gene=F, 
#                              use.log=T) # don't scale with only two columns!
# p
# 
# p <- p+scale_fill_viridis_c(limits=c(0, 2.2), oob=scales::squish)#,
#                             breaks=c(0,1,2,3,4),
#                             labels=c('0','1','2','3','> 4'))
# p
# ggsave(paste0(plotDir, 'slide_10-1.pdf'))

# split because no good scale for all these genes
high.gene.vec <- c("BICC1", "DCDC2")
low.gene.vec <- c("FGF13", "FMNL2", "FGFR2", "SERPINE1", 'KLF6', "BCL2", "CDH6", "CREB5", "KRT23", "SOX4")

p.l <- make_expression_heatmap(seu, low.gene.vec, 'is.clust.22',
                             scale.by.gene=F, 
                             use.log=F) # don't scale with only two columns!
p.l <- p.l + scale_fill_viridis_c(limits=c(0, 5), oob=scales::squish,
                            breaks=c(0,1,2,3,4,5),
                            labels=c('0','1','2','3','4', '> 5'))+
  theme(legend.position='right', 
        plot.margin=margin(t=3, r=10, b=10, l=10)) +
  scale_x_discrete(limits=rev)
p.l
p.h <- make_expression_heatmap(seu, high.gene.vec, 'is.clust.22',
                             scale.by.gene=F, 
                             use.log=F) # don't scale with only two columns!
p.h <- p.h + 
  theme(legend.position='right', 
        axis.text.x=element_blank(), 
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(t=20, r=10, b=3, l=10)) +
  scale_x_discrete(limits=rev)

p.h
p.l

p <- plot_grid(p.h, p.l, ncol=1, align='v', rel_heights=c(0.315, 1))
p
ggsave(paste0(plotDir, 'slide_10-1+SOX4+KLF6.pdf'), 
       width=4, height=6.3)



# slide 11 ---------------------------------------------------------------------
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'))
m <- seu@meta.data

p <- DimPlot(seu, reduction='umap_harmony_t.0.1', group.by='cell.annotation')+
  theme(legend.position = 'bottom')+
  xlab('')+
  ylab('')
ggsave(paste0(plotDir, 'slide_11-1.pdf'), width=6, height=7)

p <- DimPlot(seu, reduction='umap_harmony_t.0.1', group.by='Disease.status')+
  theme(legend.position = 'bottom', 
        legend.text=element_text(size=6))+
  xlab('')+
  ylab('')
ggsave(paste0(plotDir, 'slide_11-2.pdf'), width=6, height=7)


# slide 12 ---------------------------------------------------------------------
hep.seu <- readRDS(file=paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0_c5_subclusters.rds'))

plot.endstage=T
# switch which seu obj to load whether want end stage only, or all chol-hep
if (plot.endstage) {
  chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.rds_end_stage.rds'))
  reduction='umap_harmony_t.0' 
  out.suffix='_endStage.pdf'
} else {
  chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'))
  reduction='umap_harmony_t.0.1'
  out.suffix='.pdf'
}


# add subcluster info ---
chol.hep.seu@meta.data$subcluster <- chol.hep.seu@meta.data$cell.annotation
# add hep subcluster info
subclust='0'
for (subclust in c('2', '3')) {
  barcodes <- rownames(hep.seu@meta.data[hep.seu@meta.data$clust.9.subcluster==subclust,])
  chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- 'chol. like heps'
}
unique(chol.hep.seu@meta.data$subcluster)
# add chol 5 subcluster info
for (subclust in unique(chol.seu@meta.data$cluster.5.subcluster)) {
  if (subclust=='NO') {
    next
  } else {
    barcodes <- rownames(chol.seu@meta.data[chol.seu@meta.data$cluster.5.subcluster==subclust,])
    chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- paste0('chol-5')
  }
}
# add chol subcluster 22 info as well ---
barcodes <- rownames(chol.seu@meta.data[chol.seu@meta.data$SCT_snn_harmony_t.0.1.6=="22",])
chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- paste0('chol-22')
unique(chol.hep.seu@meta.data$subcluster)


# make the UMAPs
p <- DimPlot(chol.hep.seu, reduction=reduction, group.by='subcluster',
             pt.size=0.8, order=c('chol. like heps')) +
  scale_color_manual(values=c('chol-22'='grey87',
                              'chol-5'='grey87',
                              'chol. like heps'='blueviolet',
                              'Cholangiocytes'='grey87',
                              'Hepatocytes'='grey87'))+
  theme(legend.position='bottom',
        legend.text=element_text(size=6))
p
ggsave(paste0(plotDir, 'slide_12-1', out.suffix), width=7, height=7)

# p <- DimPlot(chol.hep.seu, reduction='umap_harmony_t.0.1', group.by='subcluster') +
#   scale_color_manual(values=c('chol-22'='#00BFC4',
#                               'chol-5'='#F8766D',
#                               'chol. like heps'='grey87',
#                               'Cholangiocytes'='grey87',
#                               'Hepatocytes'='grey87'))+
#   theme(legend.position='bottom', 
#         legend.text=element_text(size=6))
# p
# ggsave(paste0(plotDir, 'slide_12-2.pdf'), width=5, height=6)
# 
# p <- DimPlot(chol.hep.seu, reduction='umap_harmony_t.0.1', group.by='subcluster') +
#   scale_color_manual(values=c('chol-22'='#00BFC4',
#                               'chol-5'='#F8766D',
#                               'chol. like heps'='#7CAE00',
#                               'Cholangiocytes'='grey87',
#                               'Hepatocytes'='grey87'))+
#   theme(legend.position='bottom', 
#         legend.text=element_text(size=6))
# p
# ggsave(paste0(plotDir, 'slide_12-3.pdf'), width=5, height=6)

p <- DimPlot(chol.hep.seu, reduction=reduction, group.by='subcluster',
             pt.size=0.8, order=c('chol-22', 'chol. like heps')) +
  scale_color_manual(values=c('chol-22'='tomato',
                              'chol-5'='grey87',
                              'chol. like heps'='blueviolet',
                              'Cholangiocytes'='grey87',
                              'Hepatocytes'='grey87'))+
  theme(legend.position='bottom', 
        legend.text=element_text(size=6))
p
ggsave(paste0(plotDir, 'slide_12-4', out.suffix), width=7, height=7)



# Slide 12b - is slide 12-like but different clusters shown ------------------------
hep.seu <- readRDS(file=paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))
chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'))

table(chol.seu@meta.data[, c('SCT_snn_harmony_t.0.0.4', 'SCT_snn_harmony_t.0.1.6')])

# add subcluster info ---
chol.hep.seu@meta.data$subcluster <- chol.hep.seu@meta.data$cell.annotation
# add hep subcluster info
clust.col='SCT_snn_harmony_t.0.0.4'
for (subclust in c('9', '7')) {
  barcodes <- rownames(hep.seu@meta.data[hep.seu@meta.data[[clust.col]]==subclust,])
  chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- paste0('hep-', subclust)
}
unique(chol.hep.seu@meta.data$subcluster)
# add chol 5 subcluster info
clust.col <- 'SCT_snn_harmony_t.0.0.4'
for (subclust in c('5', '9', '1')) {
  barcodes <- rownames(chol.seu@meta.data[chol.seu@meta.data[[clust.col]]==subclust,])
  chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- paste0('chol-', subclust)
}
unique(chol.hep.seu@meta.data$subcluster)

p <- DimPlot(chol.hep.seu, reduction='umap_harmony_t.0.1', group.by='subcluster',
             pt.size=0.8, order=c('chol-9', 'chol-5', 'chol-1', 'hep-9', 'hep-7')) +
  scale_color_manual(values=c('Cholangiocytes'='grey87',
                              'Hepatocytes'='grey87',
                              'hep-7'='dodgerblue3',
                              'hep-9'='dodgerblue1',
                              'chol-1'='tomato',
                              'chol-5'='tomato3',
                              'chol-9'='tomato4'))+
  theme(legend.position='bottom', 
        legend.text=element_text(size=6))
p
ggsave(paste0(plotDir, 'slide_12b.pdf'), width=5, height=6)

# make same plot, for each cluster separately
col.code <- list('Cholangiocytes'='grey87',
  'Hepatocytes'='grey87',
  'hep-7'='dodgerblue3',
  'hep-9'='dodgerblue1',
  'chol-1'='tomato',
  'chol-5'='tomato3',
  'chol-9'='tomato4')

curr.clust <- 'chol-9'
pdf(paste0(plotDir, 'slide_12b_split.pdf'), w=5, h=6)
for(curr.clust in names(col.code)) {
  if (curr.clust %in% c('Hepatocytes', 'Cholangiocytes')) {
    next
  }
  curr.colors <- c('Cholangiocytes'='grey87',
                  'Hepatocytes'='grey87',
                  'hep-7'='grey87',
                  'hep-9'='grey87',
                  'chol-1'='grey87',
                  'chol-5'='grey87',
                  'chol-9'='grey87')
  curr.colors[curr.clust] <- col.code[[curr.clust]]

  p <- DimPlot(chol.hep.seu, reduction='umap_harmony_t.0.1', group.by='subcluster',
               pt.size=0.8, order=c(curr.clust)) +
    scale_color_manual(values=curr.colors)+
    theme(legend.position='bottom', 
          legend.text=element_text(size=6))
  print(p)
}
dev.off()


# Slide 12c - is slide 12-like but different clusters shown ------------------------
hep.seu <- readRDS(file=paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0.rds'))

plot.endstage=T
# switch which seu obj to load whether want end stage only, or all chol-hep
if (plot.endstage) {
  chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.rds_end_stage.rds'))
  reduction='umap_harmony_t.0' 
  out.suffix='_endStage.pdf'
} else {
  chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'))
  reduction='umap_harmony_t.0.1'
  out.suffix='.pdf'
}
table(chol.seu@meta.data[, c('SCT_snn_harmony_t.0.0.4', 'SCT_snn_harmony_t.0.1.6')])

# add subcluster info ---
chol.hep.seu@meta.data$subcluster <- chol.hep.seu@meta.data$cell.annotation
# add chol 5 subcluster info
clust.col <- 'SCT_snn_harmony_t.0.1.6'
for (subclust in c('3', '7', '22')) {
  barcodes <- rownames(chol.seu@meta.data[chol.seu@meta.data[[clust.col]]==subclust,])
  chol.hep.seu@meta.data$subcluster[rownames(chol.hep.seu@meta.data) %in% barcodes] <- paste0('chol-', subclust)
}
unique(chol.hep.seu@meta.data$subcluster)


p <- DimPlot(chol.hep.seu, reduction=reduction, group.by='subcluster',
             pt.size=0.8, order=c('chol-22', 'chol-7', 'chol-3')) +
  scale_color_manual(values=c('Cholangiocytes'='grey87',
                              'Hepatocytes'='grey87',
                              'chol-3'='tomato',
                              'chol-7'='tomato3',
                              'chol-22'='tomato4'))+
  theme(legend.position='bottom', 
        legend.text=element_text(size=6))
p
ggsave(paste0(plotDir, 'slide_12c', out.suffix), width=5, height=6)

# make same plot, for each cluster separately
col.code <- list('Cholangiocytes'='grey87',
                 'Hepatocytes'='grey87',
                 'chol-3'='tomato',
                 'chol-7'='tomato3',
                 'chol-22'='tomato4')

curr.clust <- 'chol-9'
pdf(paste0(plotDir, 'slide_12c_split', out.suffix), w=5, h=6)
for(curr.clust in names(col.code)) {
  if (curr.clust %in% c('Hepatocytes', 'Cholangiocytes')) {
    next
  }
  curr.colors <- c('Cholangiocytes'='grey87',
                   'Hepatocytes'='grey87',
                   'chol-3'='grey87',
                   'chol-7'='grey87',
                   'chol-22'='grey87')
  curr.colors[curr.clust] <- col.code[[curr.clust]]
  
  p <- DimPlot(chol.hep.seu, reduction=reduction, group.by='subcluster',
               pt.size=0.8, order=c(curr.clust)) +
    scale_color_manual(values=curr.colors)+
    theme(legend.position='bottom', 
          legend.text=element_text(size=6))
  print(p)
}
dev.off()




# Slide 13 ---------------------------------------------------------------------
seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.rds_end_stage.rds'))

p <- DimPlot(seu, reduction='umap_harmony_t.0', group.by='cell.annotation') +
  theme(legend.position='bottom')
ggsave(paste0(plotDir, 'slide_13.pdf'), width=7, height=7)


# Slide 14 ---------------------------------------------------------------------
# made in "pseudotime_based_analysis.R"

# Slide 16 ---------------------------------------------------------------------
# have already made the RNA velocity with chol-22 only. (chol-22-Velocity-Plot-cell.annotation.pdf)

# make the rds objects for the RNA velocity with combined 9.2 & 9.3 in the same format as did previously
# i.e. with the cluster info in the cell.annotation column
hep.seu <- readRDS(file=paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'))
chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'))

hep.barcodes <- rownames(hep.seu@meta.data[hep.seu@meta.data$clust.9.subcluster %in% c('2', '3'),])

chol.hep.seu@meta.data$cell.annotation[rownames(chol.hep.seu@meta.data) %in% hep.barcodes] <- 'chol.like.heps'
unique(chol.hep.seu@meta.data$cell.annotation)
saveRDS(chol.hep.seu, 
        file=paste0(projDir, seuDir, 'cluster_of_interest_objects/Aggr_Jan2022_Chol-Hep_harmony_th=0.1_chol.like.heps.rds'))

# then use this object with 7.RNA_velcocity/Scripts/run_3.Compute-RNA-velocity.R
# and run_4.Velocity.estimate-to-Velocity-plot.R


# Slide  17 --------------------------------------------------------------------
chol.hep.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds'))
hep.seu <- readRDS(file=paste0(projDir, seuDir, 'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'))
chol.seu <- readRDS(paste0(projDir, seuDir, 'Aggr_Jan2022_Chol_harmony_th=0_c5_subclusters.rds'))

################################################################################################
THIS IS NOW MIGRATED TO do_clust_o_I_GSEA_and_barplot.R until Chris+Vas decide which to use!!
########################################################################################

# # DE analysis of the "bridge" cells - (hep9.2, 9.3, & chol 22) vs 
# # chol-hep ---
# clust.of.interest.cells <- list()
# # which cells are in the chol clusters of interest
# clust.of.interest.cells[['chol-22']] <- rownames(chol.seu@meta.data)[chol.seu@meta.data$SCT_snn_harmony_t.0.1.6==22]
# # which cells are in the hep subclusters of interest
# clust.of.interest.cells[['hep-9.2']] <- rownames(hep.seu@meta.data)[hep.seu@meta.data$clust.9.subcluster=='2']
# clust.of.interest.cells[['hep-9.3']] <- rownames(hep.seu@meta.data)[hep.seu@meta.data$clust.9.subcluster=='3']
# all.barcodes.of.interest <- unlist(clust.of.interest.cells)
# 
# chol.hep.seu@meta.data$is.bridge.cluster <- chol.hep.seu@meta.data$cell.annotation
# chol.hep.seu@meta.data$is.bridge.cluster[rownames(chol.hep.seu@meta.data) %in% all.barcodes.of.interest] <-'bridge'
# chol.hep.seu@meta.data$is.bridge.cluster[!(rownames(chol.hep.seu@meta.data) %in% all.barcodes.of.interest)] <-'not.bridge'
# unique(chol.hep.seu@meta.data$is.bridge.cluster)
# 
# # sanity check that got correct cells
# DimPlot(chol.hep.seu, reduction='umap_harmony_t.0.1', group.by='is.bridge.cluster')
# 
# 
# # DGE analysis for cells in subclusters which are bridge
# 
# Idents(chol.hep.seu) <- chol.hep.seu@meta.data$is.bridge.cluster
# # DE parameters:
# logFC.th=0.5
# p.val.adj.th=0.1
# 
# curr.markers <- FindMarkers(
#   chol.hep.seu,
#   ident.1='bridge',
#   ident.2='not.bridge',
#   logfc.threshold=logFC.th,
#   min.pct=0.1,
#   base=2)
# curr.markers$gene <- rownames(curr.markers)
# all.markers <- curr.markers
# 
# # filter for markers expressed higher in the cluster
# all.markers.filt <- all.markers[all.markers[['avg_log2FC']] > 0,]
# # filter for alpha value
# all.markers.filt <- all.markers.filt[all.markers.filt[['p_val_adj']]<= p.val.adj.th,]
# 
# # all.markers.filt <- all.markers.filt[all.markers.filt[['cluster']]=="TRUE",]
# # save the Wilcox markers recovered:
# fwrite(all.markers.filt, 
#        file=paste0(plotDir, 'slide_17_DE_genes.csv'))
# 
# # do the GSEA - for the all markers
# background.genes <- get_SCT_genes(chol.hep.seu)
# curr.DE.genes <- all.markers.filt$gene
# gprofiler_results = gprofiler2::gost(curr.DE.genes,
#                                      organism='hsapiens',
#                                      custom_bg=background.genes,
#                                      sources=c('GO:BP', 'GO:MF','GO:CC','KEGG','REAC','TF','MIRNA'),
#                                      correction_method = 'fdr',
#                                      evcodes=T)
#   
# # save full results table - nb only returns significant terms
# fwrite(gprofiler_results$result,
#        file=paste0(plotDir, 'slide_17_GSEA_results.csv'))
# gsea.results <- fread(paste0(plotDir, 'slide_17_GSEA_results.csv'))
# 
# # get the terms told to plot
# gsea.o.I <- fread('/home/USSR/awc30/liver_project/DEG_analysis_clusters/GSEA_of_interest_v3.csv')
# 
# GSEA_oI <- merge(gsea.results, gsea.o.I, by.x='term_id', by.y='term.id',
#                      all.y=T)
# 
# GSEA_oI <- GSEA_oI[order(p_value, decreasing=T)]
# GSEA_oI$term.name <- factor(GSEA_oI$term.name, levels=unique(GSEA_oI$term.name))
# 
# p <- ggplot(GSEA_oI, aes(x=term.name, y=-log(p_value)))+
#   geom_bar(stat='identity', position='dodge')+
#   facet_wrap(~term.type, scales='free')+
#   #scale_x_discrete(position='top')+
#   theme_bw()+
#   xlab('')+
#   ylab('-log(p value)')+
#   coord_flip()
#   #theme(axis.text.x=element_text(angle=30, hjust=0 ))
# p
# ggsave(paste0(plotDir, 'slide_17_barplot.pdf'), width=10, height=2)
# 




rm(list=ls())

library(Seurat)
library(data.table)
library(Matrix) # required to transpose sparse matrix
library(preprocessCore) # for normalize.quantiles
source('/home/USSR/awc30/liver_project/zonation_markers/scripts/show_zonation_breakdown_functions.R')


# the seurat object to use - use the Hep only data
seu.path <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/6.Aggr-Jan2022/Aggr_Jan2022_Hep_harmony_th=0.rds'

marker.tab.path <- '/home/USSR/awc30/liver_project/zonation_markers/metadata/Zonation_markers_cleaned+HAL.csv'

plot.dir <- '/home/USSR/awc30/liver_project/zonation_markers/plots/'

# Load data -----
seu <- readRDS(seu.path)
marker.df <- fread(marker.tab.path)


# clean up marker df
marker.df.m <- melt(marker.df, measure.vars=c('Pericentral', 'Periportal'))
names(marker.df.m) <- c('zone', 'gene')
marker.df.m <- unique(marker.df.m)
marker.df.m <- marker.df.m[gene!='',]
marker.df.m <- marker.df.m[order(marker.df.m$zone),]


# CORRELATION IDEA --------------------
cor.method <- 'pearson' # 'spearman' 'pearson'

ds <- 'Healthy control'
all.Cs <- list()
for (ds in unique(seu@meta.data$Disease.status)) {
  print(paste0(ds, '...'))
  
  # Cut to data of interest (i.e. the marker genes) - use SCT normalised counts
  Idents(seu) <- seu@meta.data$Disease.status
  curr.seu <- subset(seu, idents=ds)
  
  # get expression
  expr <- curr.seu@assays$SCT@counts
  marker.expr <- get_marker_expression(expr, marker.df.m)
  
  # filter out cells which don't express any of the markers
  keep <- apply(marker.expr, 2, sum) > 0
  marker.expr.filt <- marker.expr[,keep]
  
  C <- calc_gene_correlation(marker.expr.filt, cor.method)
  C$Disease.status <- ds
  all.Cs[[ds]] <- C

}
Cs <- do.call('rbind', all.Cs)

# Plot the correlation
p <- plot_correlation_matrices(Cs, marker.df)
p
ggsave(filename=paste0(plot.dir, 'marker_correlation_plot_', cor.method, '.pdf'),
       width=9, height=6)



# SCORE BASED IDEA ---------------------

curr.seu <- subset(seu)
quantile.normalise=T

# get expression
expr <- curr.seu@assays$SCT@counts # was @counts # has been normalised for library size
marker.expr <- get_marker_expression(expr, marker.df.m)

# Qunatile normalise the expression of each of the marker genes, do the
# distribution of their expression across ALL the cells (including all disease
# stages) is the same.
if (quantile.normalise) {
  marker.expr.q <- t(normalize.quantiles(t(data.matrix(marker.expr))))
  rownames(marker.expr.q) <- rownames(marker.expr)
  colnames(marker.expr.q) <- colnames(marker.expr)
  marker.expr <- marker.expr.q
}

normalize.quantiles


# Calc. "score" based on normalised expression of markers for each zone separately
# is just the mean quantile normalised expression of each set of markers
PC.markers <- marker.df.m[zone=='Pericentral',]
PP.markers <- marker.df.m[zone=='Periportal',]
pc.marker.expr <- get_marker_expression(marker.expr, PC.markers)
pp.marker.expr <- get_marker_expression(marker.expr, PP.markers)
# pc.score <- calc_binarized_score(pc.marker.expr, PC.markers)
# pp.score <- calc_binarized_score(pp.marker.expr, PP.markers)
pc.score <- calc_sum_score(pc.marker.expr)
pp.score <- calc_sum_score(pp.marker.expr)

all.scores <- data.table(Disease.status=curr.seu@meta.data$Disease.status,
                         'pc.score'=pc.score,
                         'pp.score'=pp.score, 
                         'score.softmax'=calc.softmax(pc.score, 
                                                      pp.score))

# plot the score based results -----
pointSize=0.1

dis.order <- c('Healthy control', 'NAFLD',
               'NASH w/o cirrhosis', 'NASH with cirrhosis',
               'end stage')
all.scores$Disease.status <- factor(all.scores$Disease.status, levels=dis.order)

# plot softmax distribution:
sm.dist.p <- ggplot(all.scores, aes(x=score.softmax))+
  geom_histogram(bins=30)+
  facet_wrap(~Disease.status, scales='free_y', ncol=1)+
  xlab('softmax(PC : PP)')+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
           panel.grid.minor=element_blank())
sm.dist.p

# pc score against pp score:
max.lim <- max(max(all.scores$pp.score), max(all.scores$pc.score))
pc.vs.pp <- ggplot(all.scores, (aes(y=pp.score, x=pc.score)))+
  geom_point(size=pointSize)+
  facet_wrap(~Disease.status, ncol=1)+
  xlim(0, max.lim)+
  ylim(0, max.lim)+
  xlab('PC score')+
  ylab('PP score')+
  theme_bw()
pc.vs.pp

# each of the scores against the softmax score
pc.vs.sm <- ggplot(all.scores, aes(x=score.softmax, y=pc.score))+
  geom_point(size=pointSize)+
  facet_wrap(~Disease.status, ncol=1)+
  xlab('softmax(PC : PP)')+
  ylab('PC score')+
  theme_bw()

pp.vs.sm <- ggplot(all.scores, aes(x=score.softmax, y=pp.score))+
  geom_point(size=pointSize)+
  facet_wrap(~Disease.status, ncol=1)+
  xlab('softmax(PC : PP)')+
  ylab('PP score')+
  theme_bw()


p.grid <- plot_grid(sm.dist.p, pc.vs.pp, pc.vs.sm, pp.vs.sm, 
                    rel_widths = c(0.6, 0.6, 0.6, 0.6, 0.6),
                    nrow=1)
p.grid
ggsave(filename=paste0(plot.dir, 'marker_scores_plot.pdf'),
       width=10, height=10)


# Plot umaps of cells, colored by scores
seu@meta.data$score.softmax <- all.scores$score.softmax
seu@meta.data$log.score.pp <- log(all.scores$pp.score)
seu@meta.data$log.score.pc <- log(all.scores$pc.score)


p.sm <- FeaturePlot(seu, reduction='umap_harmony_t.0', features='score.softmax')+
  scale_color_viridis_c()
p.pp <- FeaturePlot(seu, reduction='umap_harmony_t.0', features='log.score.pp')+
  scale_color_viridis_c()
p.pc <- FeaturePlot(seu, reduction='umap_harmony_t.0', features='log.score.pc')+
  scale_color_viridis_c()

p.grid <- plot_grid(p.sm, p.pp, p.pc, nrow=1)
ggsave(filename=paste0(plot.dir, 'marker_scores_UMAP.pdf'),
       width=15, height=5)


# plot the expression of the marker genes against the "position" - the softmax score
# chris will decide whether this is of interest & if so, tidy up, and do stats
# (estimate mean and CIs)
# me <- data.table(t(marker.expr))
# all.scores <- cbind(all.scores, me)
# 
# mg <- m.genes[2]
# m.genes <- names(me)
# for (mg in m.genes) {
#   p <- ggplot(all.scores, aes_string(x='score.softmax', y=mg))+
#     geom_point(size=0.1)+
#     facet_wrap(~Disease.status, nrow=1)+
#     theme_bw()
#   p
#   ggsave(filename=paste0(plot.dir, 'gene_expression_test.png'),
#          width=15, height=5)
# }
# 
# cut()


make_gene_expression_umap <- function(seu, gene, is.log=F) {
  if (!(is.log)) {
    p <- FeaturePlot(seu, slot='counts', features=gene, reduction='umap_harmony_t.0',
                   pt.size=1, order=T)
  } else {
    p <- FeaturePlot(seu, slot='data', features=gene, reduction='umap_harmony_t.0',
                     pt.size=1, order=T)
  }
  p <- p+
    xlab('UMAP 1')+
    ylab('UMAP 2')
  
  if (is.log) {
    p <- p+scale_color_gradient(low='lightyellow1', 'high'='blue', 
                         na.value='grey87',
                         limits=c(0.1, max(seu@assays$SCT@data[gene,])))
    p <- p+labs(color='log expression')
  } else {
    p <- p+scale_color_gradient(low='lightyellow1', 'high'='blue', 
                                na.value='grey87',
                                limits=c(0.1, max(seu@assays$SCT@counts[gene,])))
    p <- p+labs(color='expression')
  }
  #p
  return(p)
}


plot_violinPlot <- function(chol.hep.seu, ordered.markers, ordered.subclusters) {
  
  chol.hep.seu@meta.data$barcodes <- rownames(chol.hep.seu@meta.data)
  meta <- chol.hep.seu@meta.data[, c('barcodes', 'subcluster')]
  
  C <- GetAssayData(object=chol.hep.seu, assay='SCT', slot='data') # 'data' slot is log transformed counts
  
  C <- C[ordered.markers, ]
  C.df <- data.frame(t(C))
  C.df$barcodes <- rownames(C.df)
  C.df <- merge(meta, C.df, by='barcodes')
  C.df <- data.table(C.df)
  C.df.m <- melt(C.df, id.vars=c('barcodes', 'subcluster'), 
                 variable.name='gene', value.name='logCounts')
  C.df.m$gene <- factor(C.df.m$gene, levels=ordered.markers)
  C.df.m$subcluster <- factor(C.df.m$subcluster, 
                              levels=ordered.subclusters)
  
  p <- ggplot(C.df.m, aes(x=gene, y=logCounts, fill=subcluster))+
    #geom_boxplot()+
    geom_violin(draw_quantiles = c(0.5), scale='width', kernel="gaussian")+
    xlab('')+
    ylab('log(counts)')+
    theme_bw()
  p
  return(p)
}

# seu.filt <- seu
# umap.reduction <- 'umap_harmony_t.0'
# group.by.col <- 'Disease.status'
# man.cols <- NULL

make_UMAP_plot <- function(seu.filt,
                           umap.reduction,
                           group.by.col,
                           man.cols=NULL) {
  p <- DimPlot(seu.filt, 
               reduction=umap.reduction, 
               group.by=group.by.col)
  
  p <- p+ggtitle('') +
    xlab('UMAP1')+
    ylab('UMAP2')+
    #scale_colour_hue(direction=-1)+
    theme(axis.text=element_blank(),
          axis.ticks = element_blank())
  
  if (!(is.null(man.cols))) {
    p <- p+scale_color_manual(values=man.cols)
  }
  
  #p
  p.noleg <- p+theme(legend.position = 'none')
  #p.noleg
  return(list(p, p.noleg))
}


make_facet_umap <- function(seu, umap.reduction, 
                            group.by.col, 
                            facet.col, 
                            nrow=2) {
  Idents(seu) <- seu@meta.data[[facet.col]]
  
  facets <- unique(seu.filt@meta.data[[facet.col]])
  f <- levels(facets)[1]
  leg.list <- list()
  noleg.list <- list()
  for (f in levels(facets)) {
    curr.seu <- subset(seu, idents=f)
    L <- make_UMAP_plot(curr.seu, 
                        umap.reduction,
                        group.by.col)
    
    leg.list[[f]] <- L[[1]]+ggtitle(f)
    noleg.list[[f]] <- L[[2]]+ggtitle(f)
  }
  leg.plot <- plot_grid(plotlist=leg.list, nrow=nrow)
  noleg.plot <- plot_grid(plotlist=noleg.list, nrow=nrow)
  
  return(list(leg.plot, noleg.plot))
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


make_expression_heatmap <- function(seu, gene.vec, x.axis.col, x.order=NULL, 
                                    scale.by.gene=F, 
                                    use.log=F) {
  # gene vec specifies genes to use and plot order
  # x.axis.col specifies metadata col to use for x axis
  # if provided, x.order is x plot order
  # scale.by.gene whether to scale() each gene's expression in plot
  
  if (use.log & scale.by.gene) {
    stop('scale.by.gene AND use.log shouldnt be true!')
  }
  
  # get gene expression of interest
  if (use.log) {
    expr <- data.frame(t(log(seu@assays$RNA@data[gene.vec, ]+1)))
  } else {
    expr <- data.frame(t(seu@assays$RNA@data[gene.vec, ]))
  }
  expr$barcode <- rownames(expr)
  
  # get metadata of interest
  m <- seu@meta.data
  m$barcode <- rownames(m)
  meta <- m[, c('barcode', x.axis.col)]
  names(meta) <- c('barcode', 'x.axis.col')
  
  # combine
  data.dt <- data.table(merge(expr, meta, by='barcode'))
  data.dt.m <- melt(data.dt, id.vars=c('barcode', 'x.axis.col'),
                    variable.name='gene', value.name='expression')
  
  # calc. mean
  plt.dt <- data.dt.m[, .('avg.expression'=mean(expression)), 
                      by=.(x.axis.col, gene)]
  
  # scale by gene if required
  if (scale.by.gene) {
    plt.dt[, scale.expr:=scale(avg.expression), by=.(gene)]
    plt.dt$avg.expression <- plt.dt$scale.expr
  }
  
  # set x, y axis order
  y.order.rev <- rev(gene.vec)
  plt.dt$gene <- factor(plt.dt$gene, levels=y.order.rev)
  if (!(is.null(x.order))) {
    plt.dt[[x.axis.col]] <- factor(plt.dt[[x.axis.col]], levels=x.order)
  }
  
  p <- ggplot(plt.dt, aes(x=x.axis.col, y=gene, fill=avg.expression))+
    geom_tile()+
    xlab('')+
    ylab('')+
    labs(fill='expression')+
    theme_classic()+
    scale_fill_viridis_c()+
    theme(legend.position='bottom')
  if(scale.by.gene) {
    p <- p+labs(fill='scaled expression')
  } else if (use.log) {
    p <- p+labs(fill='log expression')
  } else {
    p <- p+labs(fill='expression')
  }
  
  p  
}

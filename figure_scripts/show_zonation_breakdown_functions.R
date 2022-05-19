get_marker_expression <- function(expr, marker.df.m) {
  markers <- marker.df.m$gene
  marker.expr <- expr[rownames(expr) %in% markers, ]
  dim(marker.expr)
  return(marker.expr)
}

calc.softmax <- function(v1, v2) {
  e.v1 <- exp(v1)
  e.v2 <- exp(v2)
  sum.e <- e.v1+e.v2
  v1.prob <- e.v1 / sum.e
  return(v1.prob)
}



calc_sum_score <- function(pc.marker.expr) {
  # naive - score is just average expression of the markers
  s <- apply(pc.marker.expr, 2, mean)
  return(s)
}


calc_binarized_score <- function(pc.marker.expr, PC.markers) {
  # calculate the proportion of the markers for the celltype expressed > 0
  pc.marker.expr[pc.marker.expr > 0] <- 1
  pc.score <- apply(pc.marker.expr, 2, sum)
  pc.score.sc <- (1+pc.score) / (1+nrow(PC.markers))
  
  return(pc.score.sc)
}


calc_gene_correlation <- function(marker.expr, cor.method){
  
  marker.expr.t <- t(marker.expr)
  dim(marker.expr.t)
  C <- cor(as.matrix(marker.expr.t), method=cor.method) # calc correlation betweem marker genes
  
  C.dt <- data.frame(C)
  C.dt$gene <- rownames(C.dt)
  C.dt.m <- melt(data.table(C.dt), id.vars='gene')
  C.dt.m$variable <- as.character(C.dt.m$variable)
  names(C.dt.m) <- c('gene1', 'gene2', 'corr')
  
  C.dt.m$corr[C.dt.m$gene1==C.dt.m$gene2] <- NA
  return(C.dt.m)
}


plot_correlation_matrices <- function(Cs, marker.df) {
  Cs$gene1 <- factor(Cs$gene1, levels=marker.df.m$gene)
  Cs$gene2 <- factor(Cs$gene2, levels=marker.df.m$gene)
  Cs$Disease.status <- factor(Cs$Disease.status, 
                              levels=c('Healthy control', 
                                       'NAFLD', 
                                       'NASH w/o cirrhosis', 
                                       'NASH with cirrhosis',
                                       'end stage'))
  
  C.limit <- max(abs(na.omit(Cs$corr)))
  p <- ggplot(Cs, aes(x=gene1, y=gene2, fill=corr))+
    geom_tile()+
    facet_wrap(~Disease.status)+
    scale_fill_gradient2(low='red', mid='white', high='blue',
                         limits=c(-1*C.limit, C.limit))+
    #scale_fill_gradient2(low='red', mid='white', high='blue')+
    theme_bw()+
    xlab('')+
    ylab('')+
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  return(p)
}

make_QC_plots <- function(seu, group.by.col, umap.reduction) {
  
  meta.df <- seu@meta.data
  
  p1 <- ggplot(meta.df, aes_string(x=group.by.col, y='nCount_RNA'))+
    geom_violin()
  p2 <- ggplot(meta.df, aes_string(x=group.by.col, y='nFeature_RNA'))+
    geom_violin()
  #p2
  p3 <- ggplot(meta.df, aes_string(x=group.by.col, y='percent.mt.RNA'))+
    geom_violin()
  #p3
  p4 <- ggplot(meta.df, aes_string(x=group.by.col, y='percent.rp.RNA'))+
    geom_violin()
  #p4
  
  p5 <- ggplot(meta.df, aes(x=log(nCount_RNA), y=log(nFeature_RNA) ))+
    geom_point(aes_string(color=group.by.col))+
    geom_abline(slope=1)
  
  p6 <- p5+facet_wrap(as.formula(paste('~', group.by.col)))+
    theme(legend.position='none')
  #p6
  
  p6.1 <- FeaturePlot(seu, reduction=umap.reduction,
                       feature='nCount_RNA')+scale_color_viridis_c()
  
  p6.2 <- FeaturePlot(seu, reduction=umap.reduction,
                       feature='nFeature_RNA')+scale_color_viridis_c()
  
  p6.3 <- FeaturePlot(seu, reduction=umap.reduction,
                       feature='percent.mt.RNA')+scale_color_viridis_c()
  
  p6.4 <- FeaturePlot(seu, reduction=umap.reduction,
                       feature='percent.rp.RNA')+scale_color_viridis_c()
  
  p7 <- ggplot(meta.df, aes_string(x=group.by.col, y='S.Score'))+
    geom_violin()
  
  p8 <- ggplot(meta.df, aes_string(x=group.by.col, y='G2M.Score'))+
    geom_violin()
  
  p9 <- ggplot(meta.df, aes(x=S.Score, y=G2M.Score))+
    geom_point(aes_string(color=group.by.col), size=1)
    
  p10 <- p9 + facet_wrap(as.formula(paste('~', group.by.col)))+
    theme(legend.position='none')
  
  p11 <- FeaturePlot(seu, reduction=umap.reduction,
                     feature='S.Score')
  
  p12 <- FeaturePlot(seu, reduction=umap.reduction,
                     feature='G2M.Score')
  
  p13 <- DimPlot(seu, reduction=umap.reduction, 
                 group.by='Phase')
  
  p14 <- DimPlot(seu, reduction=umap.reduction, 
                 group.by='Phase', split.by='Phase', ncol=2)
  
  p15 <- DimPlot(seu, reduction=umap.reduction, 
                 group.by='Disease.status', split.by='Patient.ID', ncol=6)
  
  return(list(p1, p2, p3, p4, p5, p6, p6.1, p6.2, p6.3, p6.4,
              p7, p8, p9, p10, p11, p12, p13, p14, p15))
}


get_inSeu_paths <- function(dataset.dir, pattern) {
  # get the paths to the files in dataset.dir, which match the pattern
  
  all.files <- list.files(dataset.dir, pattern=pattern)
  if (length(all.files) == 0) {
    stop('No files in seu.dir match the seu.file.stem pattern specified')
  } else {
    cat(paste0('Input data files are: \n ', paste0(all.files, collapse='\n ')))
  }
  all.paths <- paste0(dataset.dir, all.files)
  return(all.paths)
}


get_SCT_genes <- function(seu) {
  # for big matrix, this breaks and too slow anyway. Probably just getting 
  
  # E <- seu@assays$RNA@counts
  # total.expression <- apply(E, 1, sum)
  # expressed <- total.expression[total.expression > 0]
  # 
  # return(names(expressed))
  
  return(rownames(seu@assays$SCT@counts))
}


load_data_and_clusters <- function(sp, cluster.col.pattern) {
  seu <- readRDS(sp)
  m <- seu@meta.data
  
  # nb, always a vector returned by grep
  clust.col <- grep(pattern=cluster.col.pattern, x=colnames(m), value=T)
  if (length(clust.col) != 1) {
    stop(paste0('cluster.col.pattern doesnt uniquely ID a cluster column. Identifies:\n', clust.col))
  }
  
  cat(paste0('using \n', clust.col, '\n as clusters for DEG analysis.'))
  Idents(seu) <- seu@meta.data[[clust.col]]
  return(list(seu, clust.col))
}

findAllBulkedMarkers <- function(seu, clust.col, bulk.by.col) {
  results.list <- list()
  seu@meta.data[[clust.col]] <- as.character(seu@meta.data[[clust.col]])
  curr.group <- unique(seu@meta.data[[clust.col]])[4]
  for (curr.group in unique(seu@meta.data[[clust.col]])) {
    
    curr.markers <- do_bulked_DE(seu, 
                                 bulk.by.col,
                                 clust.col, 
                                 curr.group)
    curr.markers$cluster <- curr.group
    curr.markers$gene <- rownames(curr.markers)
    
    results.list[[curr.group]] <- curr.markers
  }
  results <- do.call('rbind', results.list)
  return(results)    
}

# seu <- seu
# m <- seu@meta.data
# DE.comparison.col <- 'SCT_snn_harmony_t.0.0.4'
# bulk.by.col <- 'Patient.ID'
# DE.comparison.group <- 0
do_bulked_DE <- function(seu, bulk.by.col, 
                         DE.comparison.col, DE.comparison.group) {
  # using edgeR, following 
  # http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  # logFC is expression in DE.comparison.group, vs all the others, so when > 0,
  # expression is higher in the DE.comparison.group, when < 0, is lower in it.
  
  # split the object into cells in DE.comparison.group, and all others
  Idents(seu) <- seu@meta.data[[DE.comparison.col]]
  seu.group <- subset(seu, idents=DE.comparison.group)
  seu.other <- subset(seu, idents=DE.comparison.group, invert=T)
  
  # get summed gene counts across cells within each patient
  E.group <- get_aggregated_counts(seu.group, bulk.by.col)
  E.other <- get_aggregated_counts(seu.other, bulk.by.col)
  E <- data.frame(cbind(E.group, E.other)) 
  
  # E.tst <- E
  # E.tst$g <- rownames(E.tst)
  
  # get which aggregated col belongs to which test group (is E.group, then E.other)
  other.group.name <- paste0('not.', DE.comparison.group)
  data.groups <- c(rep(DE.comparison.group, ncol(E.group)),
                   rep(other.group.name, ncol(E.other)))
  data.groups <- factor(data.groups, levels=c(other.group.name, DE.comparison.group)) # levels of this factor decide which way around logFC is in final results.
  
  # meta <- unique(seu@meta.data[, c(bulk.by.col, DE.comparison.col)])
  # names(meta) <- c('id', 'group')
  # test.group <- DE.comparison.group
  # control.group <- unique(meta$group[meta$group != DE.comparison.group])
  # meta$group <- factor(meta$group, levels=c(control.group, test.group)) # levels of this factor decide which way around logFC is in final results.
  # rownames(meta) <- meta$id
  # meta <- meta[match(colnames(E), meta$id),] # need rows ids in same order as E cols
  # data.groups <- meta$group
  
  # make EdgeR object, and calc normalising factors
  d <- DGEList(counts=E, group=data.groups)
  keep <- filterByExpr(d)
  d <- d[keep,, keep.lib.sizes=F]
  d <- calcNormFactors(d) 
  design.mat <- model.matrix(~data.groups)
  d <- estimateDisp(d, design.mat) # nb design matrix used to estimate dispersion
  
  # perform the likelihood ratio test, as this is the best performing in 
  fit <- glmFit(d, design.mat)
  lrt <- glmLRT(fit, coef=2) # coef=2 specifies to test the difference between groups
  
  # extract results
  results <- lrt$table
  results$positive.logFC.means.high.in <- lrt$comparison
  results$p_val_adj <- p.adjust(results$PValue, method='BH')
  #results$gene.id <- rownames(results)
  # re-order cols
  #results <- results[, c('gene.id', names(results[1:(ncol(results)-1)]))]
  
  results <- results[order(results$PValue), ]
  
  return(results)
}


# aggregator.col <- 'Patient'
get_aggregated_counts <- function(seu, aggregator.col) {
  # get summed counts of each gene. Summed within each identity in aggregator.col
  
  expr <- seu@assays$RNA@counts
  
  aggregated.counts <- list()
  a <- unique(seu@meta.data[[aggregator.col]])[1]
  for (a in unique(seu@meta.data[[aggregator.col]])) {
    curr.barcodes <- rownames(seu@meta.data)[seu@meta.data[[aggregator.col]]==a]
    curr.expr <- expr[, curr.barcodes]
    aggregated <- apply(as.matrix(curr.expr), 1, sum)
    aggregated.counts[[a]] <- aggregated
  }
  # convert list to matrix - each column is aggregated counts for a patient
  ag.counts.m <- matrix(unlist(aggregated.counts), ncol=length(aggregated.counts))
  colnames(ag.counts.m) <- names(aggregated.counts)
  rownames(ag.counts.m) <- rownames(expr)
  
  return(ag.counts.m)
}


get_curr_outDir <- function(table.dir, sp) {
  curr.table.dir <- paste0(table.dir, basename(sp))
  curr.table.dir <- str_replace_all(curr.table.dir, '\\.rds', '')
  curr.table.dir <- paste0(curr.table.dir, '/')
  return(curr.table.dir)
}


compare_DE_results <- function(all.markers.filt, bulk.markers.filt, cluster.col) {
  # compare similarity of the DE results by 2 different methods. 
  # Each method has done DE using cluster.col as the groups DE was tested over.
  
  # check expected clusters in both DFs
  stopifnot('gene' %in% colnames(all.markers.filt))
  stopifnot('gene' %in% colnames(bulk.markers.filt))
  stopifnot(cluster.col %in% colnames(all.markers.filt))
  stopifnot(cluster.col %in% colnames(bulk.markers.filt))
  
  # screws if up if is factors
  all.markers.filt[[cluster.col]] <- as.character(all.markers.filt[[cluster.col]])
  bulk.markers.filt[[cluster.col]] <- as.character(bulk.markers.filt[[cluster.col]])
  
  
  all.clusters <- unique(c(all.markers.filt[[cluster.col]],
                           bulk.markers.filt[[cluster.col]]))
  
  all.dfs.list <- list()         
  curr.c <- all.clusters[1]                
  for (curr.c in all.clusters) {
    all.markers <- all.markers.filt[['gene']][all.markers.filt[[cluster.col]]==curr.c]
    bulk.markers <- bulk.markers.filt[['gene']][bulk.markers.filt[[cluster.col]]==curr.c]
    both.markers <- intersect(all.markers, bulk.markers)
    
    curr.df <- data.frame('cluster'=curr.c,
                          'wilcox'=length(all.markers),
                          'bulk'=length(bulk.markers),
                          'intersect'=length(both.markers))
    all.dfs.list[[curr.c]] <- curr.df
  }
  all.dfs <- do.call('rbind', all.dfs.list)
  all.dfs.m <- melt(data.table(all.dfs), id.vars=c('cluster'))
  
  # make summary plot
  p <- ggplot(all.dfs.m, aes(x=cluster, y=value, fill=variable))+
    geom_bar(stat='identity', position='dodge')+
    ylab('# DE genes')+
    theme_bw()
  
  return(p)
}


# test.use='wilcox'
# base=2
# return.thres=p.val.adj.th
# min.pct=0.1
# logfc.threshold=0.5
FindAllPairwiseMarkers <- function(seu, 
                                   logfc.threshold, 
                                   test.use='wilcox',
                                   min.pct,
                                   return.thres,
                                   base) {
  
  clusters <- as.numeric(unique(Idents(seu)))-1 # because subsetting on idents behaviour is ridiculous
                                                # if ci/cj are factors
  ci <- clusters[1]
  cj <- clusters[2]
  out.list <- list()
  for (ci in clusters) {
    for (cj in clusters) {
      if (ci==cj) {
        next
      } else {
        curr.seu <- subset(seu, idents=c(ci, cj))
        curr.markers <- FindMarkers(curr.seu,
                                    ident.1=ci,
                                    ident.2=cj,
                                    logfc.threshold=logfc.threshold,
                                    test.use=test.use,
                                    min.pct=min.pct,
                                    return.thres=return.thres,
                                    base=base)
        curr.markers$ident1 <- ci
        curr.markers$ident2 <- cj
        curr.markers$pos.logFC.means.higher.in <- ci
        curr.markers$gene <- rownames(curr.markers)
        out.list[[paste0(ci, '-', cj)]] <- curr.markers
      }
    }
  }
  out <- do.call('rbind', out.list)
  return(out)
}



# group.chol <- 'Disease.status'
# g1 <- 'NAFLD'
# g2 <- 'NASH w/o cirrhosis'
# cluster.col <- 'SCT_snn_harmony_t.0.0.4'
# cluster.of.interest=5
# chol.m <- hep.m.c16
# dis.1 <- 'Healthy control'
# dis.2 <- 'end stage'
# cluster.col <- 'new.name'
# cluster.of.interest <- 'c16' # "chol. like hep." # 'c16'
make_pairwise_cluster_proportion_comparison <- function(chol.m,
                                                        dis.1, dis.2,
                                                        cluster.col,
                                                        cluster.of.interest) {
  
  if (dis.1==dis.2) {
    return(NA)
  }
  
  # calculate p-value for a difference in the proportion of cells in the 
  # cluster of interest between the two disease.status' specified.
  
  # uses mixed model with PatientID as random effect.
  
  df <- chol.m[chol.m$Disease.status %in% c(dis.1, dis.2)]
  df$success <- df[[cluster.col]]==cluster.of.interest
  
  #table(df[, c('success', 'Disease.status')])
  
  
  if (length(unique(df$success)) == 1) {
    return(NA)
  } else {
    model <- glmer(success ~ Disease.status + (1|Patient.ID), 
                   data=df,
                   family=binomial,
                   control=glmerControl(optimizer='bobyqa',
                                        optCtrl=list(maxfun=2e5)))
    coeffs <- summary(model)$coefficients
    p.vals <- coeffs[, colnames(coeffs)=='Pr(>|z|)']
    p.dis <- p.vals[2]
    
    return(p.dis)
  }
}

# chol.m <- hep.m
# cluster.col <- 'SCT_snn_harmony_t.0.0.4'
# cluster.id <- 9
get_proportion_of_cells_in_cluster <- function(chol.m, cluster.col, cluster.id) {
  chol.m$is.CoI <- F
  chol.m$is.CoI[chol.m[[cluster.col]]==cluster.id] <- T
  dis.summary <- chol.m[, .('num.cells'=.N, 'num.cluster.of.interest'=sum(is.CoI)),
                        by=.(Disease.status)]
  dis.summary$proportion <- dis.summary$num.cluster.of.interest / dis.summary$num.cells
  
  return(dis.summary)
}


calc_all_pairwise_disease_pvals <- function(chol.m, cluster.col, cluster.of.interest) {
  results <- list()
  for (d1 in unique(chol.m$Disease.status)) {
    for (d2 in unique(chol.m$Disease.status)) {
      p.diff <- make_pairwise_cluster_proportion_comparison(chol.m,
                                                            d1, d2, 
                                                            cluster.col,
                                                            cluster.of.interest)
      results[[paste0(d1, '-', d2)]] <- data.frame('dis.1'=d1,
                                                   'dis.2'=d2,
                                                   'p.diff'=p.diff)
    }
  }
  p.df <- do.call('rbind', results)
  
  dis.order <- c('Healthy control', 'NAFLD', 'NASH w/o cirrhosis', 
                 'NASH with cirrhosis', 'end stage')
  p.df$dis.1 <- factor(p.df$dis.1, levels=dis.order)
  p.df$dis.2 <- factor(p.df$dis.2, levels=dis.order)
  
  return(p.df)
}


make_barplot <- function(cluster.prop.dt, titleStr) {
  dis.order <- c('Healthy control', 'NAFLD', 'NASH w/o cirrhosis', 
                 'NASH with cirrhosis', 'end stage')
  
  cluster.prop.dt$Disease.status = factor(cluster.prop.dt$Disease.status, 
                                          levels=dis.order)          
  p.bar <- ggplot(cluster.prop.dt, aes(x=Disease.status, y=proportion))+
    geom_bar(stat='identity')+
    theme_bw()+
    xlab('')+
    ylab('proportion of cells')+
    ggtitle(titleStr)+
    theme(panel.grid=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1, vjust=1))
  return(p.bar)
}


plot_probability <- function(pVal.df) {
  pVal.df$label <- signif(pVal.df$p.diff, digits=3)
  p.prob.plot <- ggplot(pVal.df, aes(x=dis.1, y=dis.2, fill=p.diff))+
    geom_tile()+
    xlab('')+
    ylab('')+
    geom_text(aes(label=formatC(label, 
                                format='e',
                                digits=2)), color='cyan',
              size=2)+
    scale_fill_viridis_c(direction=-1)+
    theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
  return(p.prob.plot)
}


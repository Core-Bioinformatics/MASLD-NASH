### Functions ### --------------------------------------------------------------
# seu <- seu.sub
# mon <- mon.sub

get_expression_corr <- function(mon, GoI, corr.method='spearman') {
  
  pt.matrix <- normalized_counts(mon,
                                 norm_method='log')[, order(pseudotime(mon))]
  
  row.o.i <- pt.matrix[GoI,]
  out.genes <- c()
  out.cs <- c()
  for (i in 1:nrow(pt.matrix)) {
    if (i %% 1 == 0) {
      print(paste0(i, '/', nrow(pt.matrix)))
    }
    curr.gene <- rownames(pt.matrix)[i]
    curr.row <- pt.matrix[i,]
    c <- cor(row.o.i, curr.row, method=corr.method)
    out.genes <- c(out.genes, curr.gene)
    out.cs <- c(out.cs, c)
  }
  out <- data.frame('gene.of.interest'=GoI,
                    'comparison.gene'=out.genes,
                    'correlation'=out.cs)
  return(out)
}

plot_celltype_pseudotime <- function(seu, mon) {
  
  #traj.coord<- data.frame(mon@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])
  traj.coord <- data.frame(pseudotime(mon))
  names(traj.coord) <- 'pseudotime'
  traj.coord$barcode <- rownames(traj.coord)
  
  # combine with metadata
  m <- seu@meta.data
  m$barcode <- rownames(m)
  m <- merge(m, traj.coord, by='barcode')
  
  # order celltypes by median pseudotime
  tmp <- data.table(m)
  tmp[, med.time:=median(pseudotime), by=cell.annotation]
  tmp <- unique(tmp[, c('cell.annotation', 'med.time')])
  tmp <- tmp[order(tmp$med.time)]
  
  plot.df <- m
  plot.df <- plot.df[!(plot.df$cell.annotation=='unknown'),]
  plot.df$cell.annotation <- factor(plot.df$cell.annotation, levels=tmp$cell.annotation)
  p <- ggplot(plot.df, aes(x=cell.annotation, y=pseudotime, fill=cell.annotation))+
    geom_boxplot(outlier.size=0.5)+
    xlab('')+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
          legend.position = 'none')
  return(p)
}

# reduction='umap_harmony_t.0.1'
# mon <- mon.sub
prepare_counts_for_plotting <- function(mon, DE.genes, N=50) {
  # get counts for DE.genes, fit splines through them, scale expression, 
  # and subset to take every Nth cell
  
  print('subsetting data...')
  # if(normalise==T) {
  #   # does the built in scaling using only the subset of cells
  #   print('doing built in scaling and taking log...')
  #   # size-factor normalized and log-transformed expression matrix
  pt.matrix <- normalized_counts(mon,
                                  norm_method='log')[na.omit(match(DE.genes, rownames(rowData(mon)))),
                                                    order(pseudotime(mon))]
  
  

  # # M <- exprs(mon)
  # print('not doing built in scaling, but taking log...')
  # pt.matrix <- exprs(mon)[na.omit(match(DE.genes, rownames(rowData(mon)))), 
  #                         order(pseudotime(mon))]
  # 
  # # scale
  # if (normalise==T) {
  #   print('scaling using the gene expression across the full set of cells ...')
  #   # pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  #   full.expr <- seu@assays$SCT@counts
  #   full.expr <- full.expr[DE.genes,]
  #   
  #   for (g in rownames(full.expr)) {
  #     pt.matrix[g, ] <- do.scaling.using.all.cells(pt.matrix[g,],
  #                                                  full.expr[g,])
  #   }
  #   
  # } else {
  #   print('not scaling')
  # }
  #   
  # #print('taking log...')
  # #pt.matrix <- log(pt.matrix + 1)
  
  # save the colnames, as removed when calc splines
  cellOrder <- colnames(pt.matrix)
  # convert to spline of expression so easier to see trends + compare + subset
  print('fitting splines...')
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=15)$y}))
  
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  
  
  #dim(pt.matrix)
  colnames(pt.matrix) <- cellOrder
  # subset the data (stupid to plot all cells as columns), only keep every Nth column
  print('subsetting cells...')
  pt.matrix.sub <- pt.matrix[, seq_len(ncol(pt.matrix)) %% N == 0]
  
  return(pt.matrix.sub)
}

do.scaling.using.all.cells <- function(subV, fullV) {
  avg <- mean(fullV)
  stdev <- sd(fullV)
  
  return((subV - avg) / stdev)
}

calc_clustered_gene_order <- function(pt.matrix.sub) {
  D <- dist(pt.matrix.sub)
  hcl <- hclust(D)
  geneOrder <- hcl$labels[hcl$order]
  return(geneOrder)
}

make_meta_data_df <- function(mon, seu) {
  # get plotting metadata for each cell (pseudotime, name and celltype)
  
  pseudotime <- data.frame(pseudotime(mon))
  names(pseudotime) <- 'pseudotime'
  pseudotime$barcode <- rownames(pseudotime)
  m <- seu@meta.data
  m$barcode <- rownames(m)
  m <- merge(m, pseudotime, by='barcode')
  
  plt.metadata <- m[, c('barcode', 'cell.annotation', 'pseudotime')]
  plt.metadata$barcode <- str_replace_all(plt.metadata$barcode, '-', '.')
  
  return(plt.metadata)
}

plot_genes_of_interest <- function(mon.sub, DE.genes, font.size) {
  
  # get normalised, spline modelled, scaled, subset cells matrix of gene x cells
  pt.matrix.sub <- na.omit(prepare_counts_for_plotting(mon.sub,
                                                       DE.genes, 
                                                       N=1))
  
  # convert expression matrix to df for plotting
  plt.df <- data.frame(pt.matrix.sub)
  plt.df$gene <- rownames(plt.df)
  plt.dt <- data.table(plt.df)
  plt.dt.m <- melt(plt.dt, id.vars='gene')
  
  # cell metadata for plotting
  plt.metadata <- make_meta_data_df(mon.sub, seu.sub)
  
  # combine expression, and metadata
  plt.dt.m <- merge(plt.dt.m, plt.metadata, by.y='barcode', by.x='variable')
  font.size=font.size
  
  # order cells and genes
  # hclust to calculate the gene order for plotting 
  gene.order <- calc_clustered_gene_order(pt.matrix.sub)
  
  tmp <- unique(plt.dt.m[, c('variable', 'pseudotime')])
  tmp <- tmp[order(tmp$pseudotime),]
  plt.dt.m$variable <- factor(plt.dt.m$variable, levels=tmp$variable)
  plt.dt.m$gene <- factor(plt.dt.m$gene, levels=gene.order)
  # just metadata for same
  meta.dt.m <- unique(plt.dt.m[, c('variable', 'cell.annotation', 'pseudotime')])
  meta.dt.m$variable <- factor(meta.dt.m$variable, levels=tmp$variable)
  
  p <- make_pseudotime_expression_heatmap(plt.dt.m, meta.dt.m, font.size)
  return(p)
}


make_pseudotime_expression_heatmap <- function(plt.dt.m, meta.dt.m, gene.text.size) {
  
  p.heatmap <- ggplot(plt.dt.m, aes(x=variable, y=gene, fill=value, color=value))+
    geom_tile()+
    scale_fill_viridis_c()+
    scale_color_viridis_c()+
    xlab('')+
    ylab('')+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=gene.text.size))
  #p.heatmap
  
  # make rugplots showing pseudotime and cell type
  p.meta.1 <- ggplot(meta.dt.m, aes(x=variable, y=0, color=pseudotime))+
    geom_point(shape='|', size=5)+
    theme_classic()+
    xlab('pseudotime')+
    ylab('')+
    theme(axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          axis.line.y=element_blank(), 
          legend.position='none')
  p.meta.1
  
  p.meta.2 <- ggplot(meta.dt.m, aes(x=variable, y=0, color=cell.annotation))+
    geom_point(shape='|', size=5)+
    xlab('cell type')+
    ylab('')+
    theme_classic()+
    theme(axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          axis.line.y=element_blank(), 
          legend.position='none')
  p.meta.2
  
  p.all <- plot_grid(p.heatmap, p.meta.1, p.meta.2, ncol=1, rel_heights = c(1, 0.05, 0.05),
                     align='v', axis='lr')
  return(p.all)
}

# gene.text.size=4
make_pseudotime_expression_heatmap_markers_only <- function(plt.dt.m, meta.dt.m, 
                                                            gene.text.size,
                                                            order.by='cluster') {
  
  plt.dt.m$ylab <- paste0(plt.dt.m$broad.celltype, ' :: ', plt.dt.m$gene)
  
  if (order.by=='cellType') {
    print('ordering by cellType!')
    plt.dt.m$broad.celltype <- factor(plt.dt.m$broad.celltype, 
                                      levels=c('Hepatocytes', 
                                               'Dedifferentiated progenitor',
                                               'Stem cell', 
                                               'Bridge',
                                               'Cholangiocytes'))
    plt.dt.m <- plt.dt.m[order(plt.dt.m$broad.celltype),]
    #plt.dt.m$gene <- factor(plt.dt.m$gene, levels=unique(plt.dt.m$gene))
    plt.dt.m$ylab <- factor(plt.dt.m$ylab, levels=unique(plt.dt.m$ylab))
    
    # re-order genes manually
    ylab.order <- c("Hepatocytes :: ALB", "Hepatocytes :: ASGR1",
                    "Hepatocytes :: ASS1", "Hepatocytes :: CYP3A4", 
                    "Hepatocytes :: GHR", "Hepatocytes :: GLUL", 
                    "Hepatocytes :: TTR", 
                    "Bridge :: KLF6", "Bridge :: SERPINE1",
                    "Bridge :: FKBP5", 
                    "Bridge :: FN1", "Bridge :: GPC6", "Bridge :: IL18",
                    "Bridge :: BICC1", "Bridge :: CDH6", 
                    "Bridge :: CREB5", "Bridge :: DCDC2", "Bridge :: FGF13", 
                    "Bridge :: FGFR2",  
                    "Bridge :: HDAC9", "Bridge :: KRT23", 
                    "Bridge :: NCAM1", "Bridge :: PTCHD4", 
                    "Bridge :: SOX4", "Bridge :: SOX6", "Cholangiocytes :: AQP1", 
                    "Cholangiocytes :: CFTR", "Cholangiocytes :: CTNND2", 
                    "Cholangiocytes :: EPCAM", "Cholangiocytes :: HNF1B", 
                    "Cholangiocytes :: KRT7", "Cholangiocytes :: SOX9")  
    plt.dt.m$ylab <- factor(plt.dt.m$ylab, levels=ylab.order)
    
    
    
  } else if (order.by=='cluster') {
    print('ordering by gene expression!')
    plt.dt.m <- plt.dt.m[order(plt.dt.m$gene),]
    plt.dt.m$ylab <- factor(plt.dt.m$ylab, levels=unique(plt.dt.m$ylab))
  }

  p.heatmap <- ggplot(plt.dt.m, aes(x=variable, y=ylab, fill=value, color=value))+
    geom_tile()+
    scale_fill_viridis_c()+
    scale_color_viridis_c()+
    xlab('')+
    ylab('')+
    labs(fill='scaled expression', color='scaled expression')+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=gene.text.size))
  p.heatmap
  

  # make rugplots showing pseudotime and cell type
  p.meta.1 <- ggplot(meta.dt.m, aes(x=variable, y=0, color=pseudotime))+
    geom_point(shape='|', size=5)+
    theme_classic()+
    xlab('pseudotime')+
    ylab('')+
    theme(axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          axis.line.y=element_blank(), 
          legend.position='none')
  p.meta.1
  
  p.meta.2 <- ggplot(meta.dt.m, aes(x=variable, y=0, color=cell.annotation))+
    geom_point(shape='|', size=5)+
    xlab('cell type')+
    ylab('')+
    theme_classic()+
    theme(axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          axis.line.y=element_blank(), 
          legend.position='none')
  p.meta.2
  
  p.all <- plot_grid(p.heatmap, p.meta.1, p.meta.2, ncol=1, rel_heights = c(1, 0.1, 0.1),
                     align='v', axis='lr')
  return(p.all)
}

# 1. Convert Seurat class object to Monocle3
convert.seurat.mon3 <- function(seu, assay){
  ## Create Monocle3 object from Seurat object
  suppressPackageStartupMessages(require(monocle3))
  mon <- new_cell_data_set(expression_data = seu@assays[[assay]]@counts,
                           cell_metadata = seu@meta.data, # cluster 
                           gene_metadata = data.frame("gene_short_name"=rownames(seu@assays[[assay]]@counts),
                                                      row.names=rownames(seu@assays[[assay]]@counts)))
  return(mon)
}

# 2. Pre-process monocle object  
# nb, don't think this actually makes any difference, as reduction.method
# is UMAP, which is the pre-computed one anyway!
preprocess.mon3 <- function(mon, dim.red.method = "PCA", num.dim = 30, norm.method = "none"){
  ## Perform normalization (if stated) and dimensionality reduction 
  mon <- preprocess_cds(mon, method = dim.red.method, 
                        num_dim = num.dim, 
                        norm_method = norm.method)
  return(mon)
}

# 3. Add Seurat UMAP embedding to monocle3 object (only has PCA currently)
# this is the space the trajectory will be calculated on (so don't have to 
# worry about batch correction etc, as already done by harmony (if pass
# harmony corrected reduction)).
add.seurat.UMAP <- function(seu, mon, reduction){
  # checks
  if(!reduction %in% names(seu@reductions)) stop (paste0("Reduction ", reduction, " is not present in object Dim Red. \n Available reductions are: ", paste(names(seu@reductions), collapse=", ")))
  reducedDims(mon)$UMAP <- Embeddings(seu, reduction)[colnames(mon), ]
  mon
}

# 4. Cluster cells - using the reduction specified in the call to this Rscript, 
# which has been saved to the "UMAP" field in step 3
cluster.cells.mon3 <- function(mon, reduction.method = "UMAP", cl.method = "leiden", seu.cl.annot){
  ## Cluster cells with Leiden algorithm
  mon = cluster_cells(cds = mon, 
                      reduction_method = reduction.method,
                      k = 30, 
                      cluster_method = cl.method)
  ## Replace clusters with previous ones in Seurat object
  # mon@clusters$UMAP[["clusters"]] <- seu.cl.annot
  return(mon)
}

# 5. Learn Graph
learn.graph.mon3 <- function(mon){
  ## Compute cell-cell graph
  mon = learn_graph(cds = mon,  
                    use_partition = F, #  use partitions calculated during cluster_cells and therefore to learn disjoint graph in each partition. If F, learns single graph
                    close_loop = F, # identify potential loop structure in the data space
                    learn_graph_control = NULL, # list of control parameters to be passed to the reversed graph embedding function. Default is NULL
                    verbose = FALSE)
  return(mon)
}

# 6. Run All monocle3 pipeline
run.all.monocle3 <- function(seu_obj, cl_annot, reduction){
  ## Run the whole monocle3 pipeline
  
  # 1. Convert Seurat  object to Monocle3
  print('converting...')
  mon = convert.seurat.mon3(seu = seu_obj, assay = "SCT")
  # 2. Pre-process monocle object  
  # changed num.dim to 100 (from 30) to consistent with monocle3 tutorial. (also 30 pretty quick).
  print('preprocessing...')
  mon = preprocess.mon3(mon = mon, dim.red.method = "PCA", num.dim = 100, norm.method = "none")
  # 3. Add Seurat UMAP embedding to monocle3 object
  print('adding umap...')
  mon = add.seurat.UMAP(seu = seu_obj, mon = mon, reduction)
  # 4. Cluster cells
  print('clustering...')
  mon = cluster.cells.mon3(mon = mon, reduction.method =  "UMAP", cl.method = "leiden", seu.cl.annot = seu_obj@meta.data[[cl_annot]])
  # 4b. once have clusters (or have used celltype as clusters), should 
  # find the markers which differentiate them (these have been computed per
  # celltype, so can be read directly).
  # Then annotate your cells according to type
  
  # 5. Learn Graph
  print('learning trajectory...')
  mon = learn.graph.mon3(mon)
  
  return(mon)
}


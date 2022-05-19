rm(list=ls())
library(data.table)
library(ggplot2)


margin_spacer <- function(x) {
  # where x is the column in your dataset
  left_length <- nchar(levels(factor(x)))[1]
  if (left_length > 8) {
    return((left_length - 8) * 4)
  }
  else
    return(0)
}


# specify where the output of cellPhoneDB can be found. 
# output of this script will be writter within this dir too.
cellPhoneDBDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/9.Cell-Cell-Interaction/1.CellPhoneDB/4.By-Disease-status/'

# specify the directory to use for in/out
expt.list=list('2.CellPhoneDB-Output-Chol-c5_vs_AllCells_NASH-w-cirrhosis_stages.rds'=c('Chol-5', 'NASH-w-cirrhosis'),
               '2.CellPhoneDB-Output-Chol-c5_vs_AllCells_NASH-wo-cirrhosis_stages.rds'=c('Chol-5', 'NASH-wo-cirrhosis'),
               '2.CellPhoneDB-Output-Chol-c5_vs_AllCells_all_stages'=c('Chol-5', 'allDisStates'),
               '2.CellPhoneDB-Output-Chol-c9_vs_AllCells_NAFLD_stages.rds'=c('Chol-9', 'NAFLD'),
               '2.CellPhoneDB-Output-Chol-c9_vs_AllCells_NASH-wo-cirrhosis_stages.rds'=c('Chol-9', 'NASH-wo-cirrhosis'),
               '2.CellPhoneDB-Output-Chol-c9_vs_AllCells_all_stages'=c('Chol-9', 'allDisStates'),
               '2.CellPhoneDB-Output-Chol-c22_vs_AllCells_NASH-w-cirrhosis_stages.rds'=c('Chol-22', 'NASH-w-cirrhosis'),
               '2.CellPhoneDB-Output-Chol-c22_vs_AllCells_all_stages'=c('Chol-22', 'allDisStates'),
               '2.CellPhoneDB-Output-Hep-c9.2,9.3_vs_AllCells_endStage'=c('Hep-c9.2,c9.3', 'endStage'),
               '2.CellPhoneDB-Output-Hep-c9.2,9.3_vs_AllCells_all_stages'=c('Hep-c9.2,c9.3', 'allDisStates'),
               '2.CellPhoneDB-Output-Chol-c3_vs_AllCells_all_stages'=c('Chol-3', 'allDisStages'), 
               '2.CellPhoneDB-Output-Chol-c3_vs_AllCells_endStage.rds'=c('Chol-3', 'endStage'),
               '2.CellPhoneDB-Output-Chol-c7_vs_AllCells_all_stages'=c('Chol-7', 'allDisStages'),
               '2.CellPhoneDB-Output-Chol-c7_vs_AllCells_endStage.rds'=c('Chol-7', 'endStage'),
               '2.CellPhoneDB-Output-Hep-c16_vs_AllCells_endStage'=c('Hep-c16', 'endStage'),
               '2.CellPhoneDB-Output-Hep-c16_vs_AllCells_all_stages'=c('Hep-c16', 'allDisStates'))


experimentDir <- names(expt.list)[1]
for (experimentDir in names(expt.list)) {
  cellGroupOfInterest <- expt.list[[experimentDir]][1]
  currDisStates <- expt.list[[experimentDir]][2]
  
  # input directory and output directory
  inDir <- paste0(cellPhoneDBDir, experimentDir, '/sample/')
  outDir <- paste0(cellPhoneDBDir, '/processedResultPlots/')
  if(!(dir.exists(outDir))) {
    dir.create(outDir)
  }
  
  # read the data ----
  # the identity columns in common across the input
  id.cols <- c('id_cp_interaction', 'interacting_pair', 
               'partner_a', 'partner_b', 'gene_a', 'gene_b',
               'secreted', 'receptor_a', 'receptor_b', 
               'annotation_strategy', 'is_integrin')
  
  # significant means
  smeans <- fread(paste0(inDir, 'significant_means.txt'))
  
  smeans.m <- unique(melt(smeans, id.vars=c(id.cols, 'rank'),
                   variable.name='interacting_cells',
                   value.name='mean_expression'))
  
  # means <- fread(paste0(inDir, 'means.txt'))
  pvals <- fread(paste0(inDir, 'pvalues.txt'))
  
  pvals.m <- unique(melt(pvals, id.vars=id.cols,
                  variable.name='interacting_cells', 
                  value.name='p_value'))
  
  # combine the expression and p-vals
  results <- merge(smeans.m, pvals.m, 
                   by=c(id.cols, 'interacting_cells'), all.y=T)
  
  rm(pvals, pvals.m, smeans, smeans.m)
  
  
  # filter and process the data ----
  
  # filter for comparisons between cell group of interest and others
  results.filt <- results[grepl(cellGroupOfInterest, results$interacting_cells),]
  
  if (nrow(results.filt)==0) {
    print(paste0('no  interactions for cluster : ', cellGroupOfInterest, ' in ', currDisStates))
    print('check this is a real cluster in the data!')
    exit()
  }
  
  # apply FDR multiple correction for p-val
  results.filt$p_value_corrected <- p.adjust(results.filt$p_value, 
                                             method='BH')
  
  # filter for significant
  results.filt <- results.filt[results.filt$p_value_corrected<=0.05,]
  
  if (nrow(results.filt)==0) {
    print(paste0('no SIGNIFICANT interactions for cluster : ', cellGroupOfInterest, ' in ', currDisStates))
    next
  }
  
  # split cell type comparison to cell_a and cell_b
  results.filt$celltype_a <- tstrsplit(as.character(results.filt$interacting_cells), '\\|')[[1]]
  results.filt$celltype_b <- tstrsplit(as.character(results.filt$interacting_cells), '\\|')[[2]]
  
  CoI_is_a <- results.filt[results.filt$celltype_a==cellGroupOfInterest,]
  CoI_is_b <- results.filt[results.filt$celltype_b==cellGroupOfInterest,]
  
  
  p1 <- ggplot(CoI_is_a, aes(y=interacting_pair, x=celltype_b))+
    geom_point(aes(color=p_value_corrected, size=mean_expression))+
    facet_wrap(~celltype_a)+
    xlab('')+
    ggtitle(paste0('first element is expressed in ', cellGroupOfInterest))+
    theme_bw()+  
    theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
          axis.text.y=element_text(size=6))
  
  p1
  
  
  p2 <- ggplot(CoI_is_b, aes(y=interacting_pair, x=celltype_a))+
    geom_point(aes(color=p_value_corrected, size=mean_expression))+
    facet_wrap(~celltype_b)+
    theme_bw()+
    xlab('')+
    ggtitle(paste0('second element is expressed in ', cellGroupOfInterest))+
    theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
          axis.text.y=element_text(size=6),
          plot.margin=margin(l=0+margin_spacer(CoI_is_b$interacting_pair)))
  p2
  
  pdf(paste0(outDir, cellGroupOfInterest,'-',  currDisStates, '_filtered_significance_plots.pdf'),
      width=10, height=15)
    print(p1)
    print(p2)
  dev.off()
}



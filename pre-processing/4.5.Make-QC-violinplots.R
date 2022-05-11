#!/usr/bin/env Rscript

# make QC violin plots of seurat object
suppressPackageStartupMessages(require("optparse"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("cowplot"))


### MAIN ----

option_list = list(
  make_option(
    c("-i", "--input_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input Seurat object'
  ), 
  make_option(
    c('-o', '--out_dir'),
    action='store',
    default=NA,
    type = 'character',
    help = 'path to directory where violinplots will go'
  )
)
  
opt <- parse_args(OptionParser(option_list=option_list))

inFile <- opt$input_seurat
#inFile <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/4.Aggr-52-samples-Aug2021/Aggr_52_samples_Aug_2021_meta.rds'
outDir <- opt$out_dir
#outDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/4.Aggr-52-samples-Aug2021/QC_violinplots_pre-filtering/'

dir.create(outDir, showWarnings = F)

seu <- readRDS(inFile)


# calc current rna, and mt % - NB this is just for plotting - object never saved!
seu[["percent.mt.RNA"]] <- PercentageFeatureSet(seu, assay="RNA",  pattern = "^MT-")
seu[["percent.rp.RNA"]] <- PercentageFeatureSet(seu, assay = "RNA", pattern = "^RPS") + PercentageFeatureSet(seu, assay = "RNA", pattern = "^RPL")


# do the plotting with GGplot for more control on order and color
meta <- seu@meta.data
meta$orig.ident <- factor(meta$orig.ident, levels=sort(unique(meta$orig.ident)))

pnCount <- ggplot(meta, aes(x=orig.ident, y=log10(nCount_RNA), color=Disease.status))+
  geom_boxplot()+
  geom_violin()+
  theme_bw()+
  xlab('')+
  ylab(expression(log[10](nCount)))+
  ggtitle('nCount_RNA')+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

pnFeature <- ggplot(meta, aes(x=orig.ident, y=log10(nFeature_RNA), color=Disease.status))+
  geom_boxplot()+
  geom_violin()+
  theme_bw()+
  xlab('')+
  ylab(expression(log[10](nFeature)))+
  ggtitle('nFeature_RNA')+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
#pnFeature

pnMT <- ggplot(meta, aes(x=orig.ident, y=percent.mt.RNA, color=Disease.status))+
  geom_boxplot()+
  geom_violin()+
  theme_bw()+
  xlab('')+
  ylab('MT %')+
  ggtitle('percent Mitochondrial')+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
#pnMT

pnRP <- ggplot(meta, aes(x=orig.ident, y=percent.rp.RNA, color=Disease.status))+
  geom_boxplot()+
  geom_violin()+
  theme_bw()+
  xlab('')+
  ylab('RP %')+
  ggtitle('percent Ribosomal')+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
#pnRP

pdf(paste0(outDir, 'violin_QCs.pdf'),
    width=20, height=10)
  print(pnCount)
  print(pnFeature)
  print(pnMT)
  print(pnRP)
dev.off()


meta$passed.non.mt <- TRUE
meta$passed.non.mt[meta$nFeature_RNA < 800] <- F
meta$passed.non.mt[meta$nCount_RNA < 1000] <- F
meta$passed.non.mt[meta$nCount_RNA > 50000] <- F
meta$passed.non.mt[meta$percent.rp.RNA > 10] <- F


# plot distribution of each of these features, and cumsum distribution
p.count <- ggplot(meta, aes(x=nCount_RNA))+
  geom_histogram(bins=100)+
  geom_vline(xintercept =1000)+
  geom_vline(xintercept =50000)+
  scale_x_continuous(trans='log10')+
  ggtitle('nCount thresholds')+
  theme_bw()

p.feature <- ggplot(meta, aes(x=nFeature_RNA))+
  geom_histogram(bins=100)+
  geom_vline(xintercept = 800)+
  scale_x_continuous(trans='log10')+
  ggtitle('nFeature threshold')+
  theme_bw()

p.mt<- ggplot(meta, aes(x=percent.mt.RNA))+
  geom_histogram(bins=100)+
  geom_vline(xintercept = 10)+
  scale_x_continuous(trans='log10')+
  ggtitle('% MT threshold')+
  theme_bw()
p.mt

p.rp <- ggplot(meta, aes(x=percent.rp.RNA))+
  geom_histogram(bins=100)+
  geom_vline(xintercept = 10)+
  scale_x_continuous(trans='log10')+
  ggtitle('% RB threshold')+
  theme_bw()
p.rp

pdf(paste0(outDir, 'distributions_and_thresholds.pdf'), width=4, height=4)
    print(p.count)
    print(p.feature)
    print(p.mt)
    print(p.rp)
dev.off()

# percent failed on MT
tmp <- meta[meta$passed.non.mt,]
num.mt.failed <- sum((tmp$percent.mt.RNA > 10) & tmp$percent.mt.RNA < 50)
num.passed <- sum(tmp$percent.mt.RNA <= 10)

p.mt.detailed <- ggplot(meta, aes(x=percent.mt.RNA, fill=passed.non.mt))+
  geom_histogram(bins=100, position='identity', alpha=0.4)+
  geom_vline(xintercept = 10)+
  scale_x_continuous(trans='log10')+
  ggtitle(paste0(num.mt.failed, ' cells failed on MT% alone (', num.passed, ' cells passed all)'))+
  theme_bw()+
  theme(legend.position='top')
p.mt.detailed

pdf(paste0(outDir, 'MT_detailed.pdf'), width=5, height=4)
print(p.mt.detailed)
dev.off()



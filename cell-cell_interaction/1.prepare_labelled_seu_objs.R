rm(list=ls())
library(Seurat)

projDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'
seuDir <- '6.Aggr-Jan2022/'

outDir <- '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/9.Cell-Cell-Interaction/1.CellPhoneDB/seuratObjects/'



all.Seu <- readRDS(file=paste0(projDir,
                              seuDir,
                              'Aggr_Jan2022_AllCells_harmony_th=0.rds'))



# chol cluster 5 in NASH with cirrhosis, NASH wo cirrhosis, and all disease stages -------------
chol.Seu <- readRDS(file=paste0(projDir,
                                seuDir,
                                'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m.Seu <- chol.Seu@meta.data
table(m.Seu[, c('Disease.status', 'SCT_snn_harmony_t.0.0.4')])

cells.o.I <- rownames(m.Seu[m.Seu$SCT_snn_harmony_t.0.0.4 %in% c('5'),])

# make the column "CellPhoneComp" which will be used for the CellPhoneDB comparisons
curr.Seu <- all.Seu
curr.Seu@meta.data$CellPhoneComp <- curr.Seu@meta.data$cell.annotation
curr.Seu@meta.data$CellPhoneComp[rownames(curr.Seu@meta.data) %in% cells.o.I] <- 'Chol-5'
# sanity check that all are cholangiocytes as expected.
table(curr.Seu@meta.data[, c('cell.annotation', 'CellPhoneComp')])

# save the full seu object
out.Seu <- curr.Seu
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c5_vs_AllCells_all_stages.rds'))

# subset to the desired stage cells - NASH
out.Seu <- subset(curr.Seu, subset=Disease.status=='NASH with cirrhosis')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c5_vs_AllCells_NASH-w-cirrhosis_stages.rds.rds'))

# subset to the desired stage cells - NASH wo cirrosis
out.Seu <- subset(curr.Seu, subset=Disease.status=='NASH w/o cirrhosis')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c5_vs_AllCells_NASH-wo-cirrhosis_stages.rds.rds'))

rm(out.Seu, curr.Seu, chol.Seu, cells.o.I)



# chol cluster 9 in NAFLD, NASH wo cirrhosis, and all disease stages -------------
chol.Seu <- readRDS(file=paste0(projDir,
                                seuDir,
                                'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m.Seu <- chol.Seu@meta.data
table(m.Seu[, c('Disease.status', 'SCT_snn_harmony_t.0.0.4')])

cells.o.I <- rownames(m.Seu[m.Seu$SCT_snn_harmony_t.0.0.4 %in% c('9'),])

# make the column "CellPhoneComp" which will be used for the CellPhoneDB comparisons
curr.Seu <- all.Seu
curr.Seu@meta.data$CellPhoneComp <- curr.Seu@meta.data$cell.annotation
curr.Seu@meta.data$CellPhoneComp[rownames(curr.Seu@meta.data) %in% cells.o.I] <- 'Chol-9'
# sanity check that all are cholangiocytes as expected.
table(curr.Seu@meta.data[, c('cell.annotation', 'CellPhoneComp')])

# save the full seu object
out.Seu <- curr.Seu
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c9_vs_AllCells_all_stages.rds'))

# subset to the desired stage cells - NAFLD
out.Seu <- subset(curr.Seu, subset=Disease.status=='NAFLD')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c9_vs_AllCells_NAFLD_stages.rds.rds'))

# subset to the desired stage cells - NASH wo cirrosis
out.Seu <- subset(curr.Seu, subset=Disease.status=='NASH w/o cirrhosis')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c9_vs_AllCells_NASH-wo-cirrhosis_stages.rds.rds'))

rm(out.Seu, curr.Seu, chol.Seu, cells.o.I)


# chol cluster 22 in  NASH with cirrhosis, and all disease stages -------------
chol.Seu <- readRDS(file=paste0(projDir,
                                seuDir,
                                'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m.Seu <- chol.Seu@meta.data
table(m.Seu[, c('Disease.status', 'SCT_snn_harmony_t.0.1.6')])

cells.o.I <- rownames(m.Seu[m.Seu$SCT_snn_harmony_t.0.1.6 %in% c('22'),])

# make the column "CellPhoneComp" which will be used for the CellPhoneDB comparisons
curr.Seu <- all.Seu
curr.Seu@meta.data$CellPhoneComp <- curr.Seu@meta.data$cell.annotation
curr.Seu@meta.data$CellPhoneComp[rownames(curr.Seu@meta.data) %in% cells.o.I] <- 'Chol-22'
# sanity check that all are cholangiocytes as expected.
table(curr.Seu@meta.data[, c('cell.annotation', 'CellPhoneComp')])

# save the full seu object
out.Seu <- curr.Seu
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c22_vs_AllCells_all_stages.rds'))

# subset to the desired stage cells - NASH with cirrosis
out.Seu <- subset(curr.Seu, subset=Disease.status=='NASH with cirrhosis')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c22_vs_AllCells_NASH-w-cirrhosis_stages.rds.rds'))

rm(out.Seu, curr.Seu, chol.Seu, cells.o.I)




# chol cluster 3 in  end stage, and all disease stages -------------
chol.Seu <- readRDS(file=paste0(projDir,
                                seuDir,
                                'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m.Seu <- chol.Seu@meta.data
table(m.Seu[, c('Disease.status', 'SCT_snn_harmony_t.0.1.6')])

cells.o.I <- rownames(m.Seu[m.Seu$SCT_snn_harmony_t.0.1.6 %in% c('3'),])



# make the column "CellPhoneComp" which will be used for the CellPhoneDB comparisons
curr.Seu <- all.Seu
curr.Seu@meta.data$CellPhoneComp <- curr.Seu@meta.data$cell.annotation
curr.Seu@meta.data$CellPhoneComp[rownames(curr.Seu@meta.data) %in% cells.o.I] <- 'Chol-3'
# sanity check that all are cholangiocytes as expected.
table(curr.Seu@meta.data[, c('cell.annotation', 'CellPhoneComp')])

# save the full seu object
out.Seu <- curr.Seu
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c3_vs_AllCells_all_stages.rds'))

# subset to the desired stage cells - NASH with cirrosis
out.Seu <- subset(curr.Seu, subset=Disease.status=='end stage')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c3_vs_AllCells_endStage.rds.rds'))

rm(out.Seu, curr.Seu, chol.Seu, cells.o.I)




# chol cluster 7 in  end stage, and all disease stages -------------
chol.Seu <- readRDS(file=paste0(projDir,
                                seuDir,
                                'Aggr_Jan2022_Chol_harmony_th=0.rds'))
m.Seu <- chol.Seu@meta.data
table(m.Seu[, c('Disease.status', 'SCT_snn_harmony_t.0.1.6')])

cells.o.I <- rownames(m.Seu[m.Seu$SCT_snn_harmony_t.0.1.6 %in% c('7'),])



# make the column "CellPhoneComp" which will be used for the CellPhoneDB comparisons
curr.Seu <- all.Seu
curr.Seu@meta.data$CellPhoneComp <- curr.Seu@meta.data$cell.annotation
curr.Seu@meta.data$CellPhoneComp[rownames(curr.Seu@meta.data) %in% cells.o.I] <- 'Chol-7'
# sanity check that all are cholangiocytes as expected.
table(curr.Seu@meta.data[, c('cell.annotation', 'CellPhoneComp')])

# save the full seu object
out.Seu <- curr.Seu
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c7_vs_AllCells_all_stages.rds'))

# subset to the desired stage cells - NASH with cirrosis
out.Seu <- subset(curr.Seu, subset=Disease.status=='end stage')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Chol-c7_vs_AllCells_endStage.rds.rds'))

rm(out.Seu, curr.Seu, chol.Seu, cells.o.I)


# 2. hepatocytes cluster 9.3 & 9.2 vs all cell types, end disease stage only----

hep.Seu <- readRDS(file=paste0(projDir,
                               seuDir,
                               'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'))
m.Seu <- hep.Seu@meta.data

cells.o.I <- rownames(m.Seu[m.Seu$SCT_snn_harmony_t.0.1.6 %in% c('16'),])

# make the column "CellPhoneComp" which will be used for the CellPhoneDB comparisons
curr.Seu <- all.Seu
curr.Seu@meta.data$CellPhoneComp <- curr.Seu@meta.data$cell.annotation
curr.Seu@meta.data$CellPhoneComp[rownames(curr.Seu@meta.data) %in% cells.o.I] <- 'Hep-c16'
# sanity check that all are cholangiocytes as expected.
table(curr.Seu@meta.data[, c('cell.annotation', 'CellPhoneComp')])
table(curr.Seu@meta.data[, c('Disease.status', 'CellPhoneComp')])

# save the full seu object
out.Seu <- curr.Seu
saveRDS(out.Seu,
        file=paste0(outDir, 'Hep-c16_vs_AllCells_all_stages.rds'))

# subset to the desired stage cells
out.Seu <- subset(curr.Seu, subset=Disease.status=='end stage')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Hep-c16_vs_AllCells_endStage.rds'))

rm(out.Seu, curr.Seu, hep.Seu, cells.o.I)



# 2. hepatocytes cluster 9.3 & 9.2 vs all cell types, end disease stage only----

hep.Seu <- readRDS(file=paste0(projDir,
                               seuDir,
                               'Aggr_Jan2022_Hep_harmony_th=0_c9_subclusters.rds'))
m.Seu <- hep.Seu@meta.data

cells.o.I <- rownames(m.Seu[m.Seu$clust.9.subcluster %in% c('2', '3'),])

# make the column "CellPhoneComp" which will be used for the CellPhoneDB comparisons
curr.Seu <- all.Seu
curr.Seu@meta.data$CellPhoneComp <- curr.Seu@meta.data$cell.annotation
curr.Seu@meta.data$CellPhoneComp[rownames(curr.Seu@meta.data) %in% cells.o.I] <- 'Hep-c9.2,c9.3'
# sanity check that all are cholangiocytes as expected.
table(curr.Seu@meta.data[, c('cell.annotation', 'CellPhoneComp')])

# save the full seu object
out.Seu <- curr.Seu
saveRDS(out.Seu,
        file=paste0(outDir, 'Hep-c9.2,9.3_vs_AllCells_all_stages.rds'))

# subset to the desired stage cells
out.Seu <- subset(curr.Seu, subset=Disease.status=='end stage')
m <- out.Seu@meta.data
table(m[, c('Disease.status', 'CellPhoneComp')])
saveRDS(out.Seu,
        file=paste0(outDir, 'Hep-c9.2,9.3_vs_AllCells_endStage.rds'))

rm(out.Seu, curr.Seu, hep.Seu, cells.o.I)


#!/usr/bin/env bash

set -eu

# project vectors from two previously calculated separate velocity calculations 
# (here using different cell sets, Chol & Hep cells) onto the same UMAP.

# specify the seurat object to load this will specify where the points are plotted
# (the UMAP)
seuDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/6.Aggr-Jan2022/'
seuObj='Aggr_Jan2022_Chol-Hep_harmony_th=0.rds_end_stage.rds' # 'Aggr_Jan2022_Chol-Hep_harmony_th=0.1.rds' # 'Aggr_Jan2022_Chol-Hep_harmony_th=0.rds_end_stage.rds'
seuFile=${seuDir}${seuObj}

# specify which seurat reduction to use for plotting
reduction='umap_harmony_t.0' # 'umap_harmony_t.0.1' #'umap_harmony_t.0' 

# specify which Chol velocity object to use
cholVeloPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/7.RNA-Velocity/2.Aggr_samples/6.Aggr-Jan2022-Chol-Hep_Cholangiocytes_endStage/Aggr_samples_Jan2022-velocitity.estimates.rds'

# specify which Hep velocity object to use
hepVeloPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/7.RNA-Velocity/2.Aggr_samples/6.Aggr-Jan2022-Chol-Hep_Hepatocytes_endStage/Aggr_samples_Jan2022-velocitity.estimates.rds'

# specify where the out plot will be saved to
plotDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/7.RNA-Velocity/2.Aggr_samples/separate_arrow_plots/'
plotFile='Chol-Hep_EndStage_vel_on_Endstage_UMAP_v2.pdf'
plotPath=${plotDir}${plotFile}

# specify whether to filter cells to end stage or not
filterToEndStage='no' #'filter' or 'no'

Rscript ./make_separate_arrow_plots.R \
 ${seuFile} \
 ${reduction} \
 ${cholVeloPath} \
 ${hepVeloPath} \
 ${plotPath} \
 ${filterToEndStage}
 
#!/usr/bin/env bash

set -eu

# Calls 3.Compute-RNA-velocity.R to calculate RNA velocity usign the specified 
# celltype/cluster

# 'Hepatocytes' or 'Cholangiocytes' or the cluster name to use for calcuating velocities
cellType=$1
# 'All' or 'endStage' - filter to this disease stage
disStage=$2

# specify the seurat object to load (nb this may still be further filtered depending on
# params passed).
seuDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/6.Aggr-Jan2022/cluster_of_interest_objects/'
seuObj='Aggr_Jan2022_Chol-Hep_harmony_th=0.1_'${cellType}'.rds'

echo "params:"
echo ${cellType}
echo ${disStage}

# specify the dir where the plot and the velo.obj will be written
VeloDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/7.RNA-Velocity/2.Aggr_samples/'
exptDir='6.Aggr-Jan2022-Chol-Hep_'${cellType}'_'${disStage}'/'
outDir=${VeloDir}${exptDir}


# --barcode_file_path : ids of barcodes want to make sure are included in 
# velocity calculation.
Rscript ./3.Compute-RNA-velocity.R \
  --input_seurat ${seuDir}${seuObj} \
  --cr_paths '../Velocyto-Paths-All-Samples_Jan2022.csv' \
  --annots 'cell.annotation,Disease.Status' \
  --sample_name 'Aggr_samples_Jan2022' \
  --subset_celltype ${cellType} \
  --subset_disease_stage ${disStage} \
  --barcode_file_path '/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/7.RNA-Velocity/Aggr_samples_Jan2022-barcodes_oI.rds' \
  --output_path ${outDir}

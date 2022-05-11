#!/usr/bin/env bash
set -eu

# calls 4.Velocity.estimate-to-Velocity-Plot.R to project (previously calculated)
# velocity results onto a UMAP

#cellTypes=( 'chol-3' 'chol-7' 'chol-22' 'chol-c5' 'hep-9.2' 'hep-9.3' )
#cellTypes=( 'chol.like.heps' 'chol-22' )
cellTypes=( 'hep-7' 'chol-9' )

for cellType in "${cellTypes[@]}"
do
  echo ${cellType}
  # cellType='chol-3'

  # which seurat object to show UMAP for
  seuProjDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'
  seuDir='6.Aggr-Jan2022/cluster_of_interest_objects/'
  seuFile='Aggr_Jan2022_Chol-Hep_harmony_th=0.1_'${cellType}'.rds'
  reduction='umap_harmony_t.0.1'
  # seuDir='6.Aggr-Jan2022/'
  # seuFile='Aggr_Jan2022_Chol-Hep_harmony_th=0.rds_end_stage.rds'
  # reduction='umap_harmony_t.0' # which reduction to use for the UMAP

  # which velocity results to project onto seurat object
  projDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/7.RNA-Velocity/2.Aggr_samples/'
  veloDirStem='/6.Aggr-Jan2022-Chol-Hep_'${cellType}'_All/'
  veloFile=${projDir}${veloDirStem}Aggr_samples_Jan2022-velocitity.estimates.rds
  
  # where output will go
  plotDir=${projDir}${veloDirStem}
  
  
  Rscript ./4.Velocity.estimate-to-Velocity-Plot.R -i ${seuProjDir}${seuDir}${seuFile} \
  -v ${veloFile} \
  -c 'cell.annotation' \
  -o ${plotDir} \
  -s ${cellType} \
  -r ${reduction}
  
done
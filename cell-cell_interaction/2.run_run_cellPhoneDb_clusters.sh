#!/usr/bin/env bash

set -eu

# the names of the .rds files with the terminal .rds removed

# allInputs=( 'Chol-c5_vs_AllCells_all_stages' 'Chol-c5_vs_AllCells_NASH-w-cirrhosis_stages.rds' 
# 'Chol-c5_vs_AllCells_NASH-wo-cirrhosis_stages.rds' 'Chol-c9_vs_AllCells_all_stages' 
# 'Chol-c9_vs_AllCells_NAFLD_stages.rds' 'Chol-c9_vs_AllCells_NASH-wo-cirrhosis_stages.rds'
# 'Chol-c22_vs_AllCells_all_stages' 'Chol-c22_vs_AllCells_NASH-w-cirrhosis_stages.rds'
# 'Hep-c9.2,9.3_vs_AllCells_all_stages' 'Hep-c9.2,9.3_vs_AllCells_endStage' )

# allInputs=( 
# 'Hep-c9.2,9.3_vs_AllCells_all_stages' 'Hep-c9.2,9.3_vs_AllCells_endStage' )

# allInputs=( 'Chol-c3_vs_AllCells_all_stages' 'Chol-c3_vs_AllCells_endStage.rds' 
# 'Chol-c7_vs_AllCells_all_stages' 'Chol-c7_vs_AllCells_endStage.rds' )

allInputs=( 'Hep-c16_vs_AllCells_all_stages' 'Hep-c16_vs_AllCells_endStage' )

for i in "${allInputs[@]}"
do
  echo ${i}
  /sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/9.Cell-Cell-Interaction/1.CellPhoneDB/Scripts/run_cellPhoneDb_clusters.sh ${i}
done
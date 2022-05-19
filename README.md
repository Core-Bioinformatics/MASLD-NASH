# NAFLD-NASH
scripts associated with NAFLD-NASH project

- `alignment` - use CellRanger to align fastq files

- `pre-processing` - read count matrices into Seurat, and preprocess gene expression data. Add cell type annotations, and apply Harmony integration. Expected run order is:
  * 1.indv-samples.Rmd
  * do_preprocessing.sh
  * generate_cellType_annotation_jan2022.R
  * apply_harmony_correction.R

  the other files are called / required by these ones.

- `RNA_velocity` - calculate RNA velocity using velocyto and plot it. Two steps: 
  1. run velocity calculations (run `run_3.Compute-RNA-velocity.sh`)
  2. plot them (run `run_4.Velocity.estimate-to-Velocity_plot.sh` or `run_make_separate_arrow_plots.sh`).

- `cell-cell_interaction` - use CellPhoneDB to identify potentially interacting ligand-receptor pairs in different cell types. Expected run order is:
 * `1.prepare_labelled_seu_objs.R`
 * `2.run_run_cellPhoneDb_clusters.sh`
 * `3.Process-cellPhoneDb-output.R`
 other scripts are required / called by these.

- `Pseudotime` - infer pseudotime in cholangiocyte and hepatocyte cells using Monocle3. Identify and plot genes differentially expressed over pseudotime.

# NAFLD-NASH
scripts associated with NAFLD-NASH project

`alignment` - use CellRanger to align fastq files

`pre-processing` - read count matrices into Seurat, and preprocess gene expression data. Add cell type annotations, and apply Harmony integration.
Expected run order is:
* 1.indv-samples.Rmd
* do_preprocessing.sh
* generate_cellType_annotation_jan2022.R
* apply_harmony_correction.R

the other files are called / required by these ones.

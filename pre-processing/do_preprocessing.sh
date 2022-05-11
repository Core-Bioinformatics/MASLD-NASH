#!/usr/bin/env bash

# wrapper script to call each of Ruben's preprocessing R scripts
# having previously run
# 1.Indv-samples.Rmd)

set -eu

# specify the path to the list of all individual seurat objects which will be
# merged (made by 1.Indv-samples-Aug2021.Rmd).
inDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/1.Indv_samples/2.Lists_objects/Batch_Dec2021_v2_+id98+testSamples/'
inFile='Dec2021_v2_+id98+testSamples.rds'
inPath=${inDir}${inFile}


# specify the path to the CR aggregated csv used to make the aggregated loupe object - used to get consistent barcodes with loupe
# aggrCsvPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/2.CR_output/Aggr_Aug2021_v2/outs/aggregation.csv'
aggrCsvPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/0.data/Aggr_Dec2021_v2_+id98+testSamples.csv'
# specify the path to the sample metadata csv
# needs to have SLX runs in "SLX" column, and Seq.id (e.g. SITTD3) in "Sequencing.ID" column
# is ok if includes info on samples which have been filtered out
metaCsvPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Scripts/Dec2021_metadata_+id98+testSamples.csv'
# metaCsvPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/Disease_status_summary_63_samples_updated__25_08_21_AC_cleaned.csv'

# specify the path to the identities of the marker genes to use for cell type annotation - (this is just used for plotting them here, no
# actual annotation done)
#markerCsvPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/2.Markers-expression/2.UMAP_markers_expression/0.Expression_UMAPs_markers/Markers-Used-In-Annotation.tsv'
markerCsvPath='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/2.Markers-expression/0.Marker_genes_sets/8.marker_genes_v4.tsv'


# define output directory where processed files will be output
aggrDir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/'
nameDir='6.Aggr-Jan2022/' # directory all output will be writen to
outDir=${aggrDir}${nameDir}
outStem=Aggr_Jan2022 # stem of the files which will be written in it


echo ${inPath}
echo ${outDir}

# setup dirs
if [[ ! -e ${outDir} ]]; then
  mkdir -p ${outDir}
fi

if [[ ! -e ${outDir}/for_celltype_annotation/ ]]; then
  mkdir -p ${outDir}/for_celltype_annotation/
fi


# Merge the list of objects to a single Seurat object
echo "merging..."
Rscript ./1.Merge-Seurat-object-list.R \
-i ${inPath} \
--list TRUE \
-o ${outDir}/${outStem}.rds
echo "done!"

# Rename the cells to be consistent between loupe aggr, and seurat aggr
echo "renaming cells..."
Rscript ./7.Rename-Cell-Names-Seurat_v2.R \
-i ${outDir}/${outStem}.rds \
--input_csv ${aggrCsvPath} \
-o ${outDir}/${outStem}_renamed.rds
echo "done!"

# Add metadata
echo "adding metadata..."
Rscript ./4.Add-Different-Metadata-CHECK_v2.R \
-i ${outDir}/${outStem}_renamed.rds \
--input_metadata ${metaCsvPath} \
-o ${outDir}/${outStem}_meta.rds
echo "done!"

# make QC violin plots
echo "making QC violinplots..."
Rscript ./4.5.Make-QC-violinplots.R \
-i ${outDir}/${outStem}_meta.rds \
-o ${outDir}/QC_violinplots_pre-filtering/
echo "done!"

# filter on MT, Rb, nFeatures, nCounts
echo "filtering..."
Rscript ./5.Seurat-filtering.R \
-i ${outDir}/${outStem}_meta.rds \
--list FALSE \
-o ${outDir}/${outStem}_filtered.rds
echo "done!"

# make QC violin plots
echo "making more QC violinplots..."
Rscript ./4.5.Make-QC-violinplots.R \
-i ${outDir}/${outStem}_filtered.rds \
-o ${outDir}/QC_violinplots_post-filtering/
echo "done!"

# remove MT, Rb genes
echo "removing mt & rb genes..."
Rscript ./6.Remove.Mito.Ribo.genes.R \
-i ${outDir}/${outStem}_filtered.rds \
--list FALSE \
-o ${outDir}/${outStem}_no_MT_RB.rds
echo "done!"

# preprocessing - SCT transform, scale, run PCA, run UMAP, find clusters
echo "preprocessing..."
Rscript ./2.Seurat-pre-processing.R \
-i ${outDir}/${outStem}_no_MT_RB.rds \
--list FALSE \
-o ${outDir}/${outStem}_pre-processed.rds
echo "done!"

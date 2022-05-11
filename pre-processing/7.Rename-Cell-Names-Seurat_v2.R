#/usr/bin/env Rscript

## Rename Cells of Seurat object. This is, replace Barcode suffix with corresponding sample number
## This is to avoid the unique cell names generate per Seurat to avoid name collission (eg. BARCODE_1_1_1_1) etc
## Sample numbers defined in in https://docs.google.com/spreadsheets/d/13W17ky2nVGZqXeC81uc-R3V1x1U8HyvDulaM2u7eu5w/edit?usp=sharing 
## in column "Actual sample num"

## Change barcodes suffix according to 'rename.vec' and replace cell names
# replaces the "_1_1_1_1..." type encoding used by seurat merge with the "_4" type
# encoding used by Cellranger.
# knows the correct mapping between names because know cellranger integer is 
# the order of the samples in the "aggregate.csv" file, and have the orig.ident
# sample names in the seurat object.

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--input_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input Seurat object, or list of Seurat objects'
  ),
  make_option(
    c("-c", "--input_csv"),
    action = "store",
    default = NA,
    type='character',
    help = 'Path to aggregated cell ranger output aggregation.csv file'
  ),
  make_option(
    c("-o", "--output_seurat"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to Seurat object or object list.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(stringr))


## Functions ## 
rename.cells.seu <- function(seu, rename.vec){
  # rename the cells using the rename.vec mappint (sample name to cellranger integer.)
  
  barcodes <- colnames(seu) 
  names(barcodes) <- seu@meta.data[["orig.ident"]]
  barcodes  <- unlist(lapply(strsplit(barcodes, "-"), function(x) x[1]))
  ## Add barcodes suffix according to rename.vec
  for (sample.name in names(rename.vec)){
    barcodes[names(barcodes) == sample.name] <- paste0(barcodes[ names(barcodes) == sample.name], "-", rename.vec[[sample.name]])
  }
  ## Add renamed Barcodes to object
  seu = RenameCells(seu, new.names = unname(barcodes)) # Unname to load vector without names
  return(seu)
}

make.sample.number.vec <- function(csv_path) {
  # make named vector of 1:numSamples, named by samples. Sample order is determined by the order or rows in csv_path
  df <- fread(csv_path)
  
  # I mistyped as SLX-19440_SITT-A2 in the CR aggregation file, so need to replace it with the correct SLX-19940_SITTA2
  # to be compatible with the (correct) naming in the seurat object
  df$library_id[df$library_id=='SLX-19440_SITT-A2'] <- 'SLX-19940_SITT-A2'
  
  # reformat df library_id to be same format as in the seurat object (rename from SLX-XXXXX_SITT-XX to SLX-XXXXX-SITTXXX)
  df$seqRun <- tstrsplit(df$library_id, '_')[1]
  df$run <- tstrsplit(df$library_id, '_')[2]
  df$run <- str_replace(df$run, '-', '')
  df$newid <- paste0(df$seqRun, '-', df$run)
  
  outvec <- 1:nrow(df)
  names(outvec) <- df$newid
  
  return(outvec)
}


sample.number.vec <- make.sample.number.vec(opt$input_csv)
#sample.number.vec <- make.sample.number.vec('/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/2.CR_output/Aggr_Aug2021_v2/outs/aggregation.csv')


## 1. Read object/s
seu = readRDS(opt$input_seurat)
#seu = readRDS('/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/4.Aggr-52-samples-Aug2021/Aggr_52_samples_Aug_2021.rds')

# check that got the same ids in both
tmp <- unique(seu[['orig.ident']])
diff <- setdiff(names(sample.number.vec), unique(tmp[['orig.ident']]))
if (length(diff) != 0) {
  stop(paste0('Discrepancy between the identities of the samples in merged seurat object 
       and aggregate.csv].\n',
              'in merged Seurat object only: ',
              setdiff(unique(tmp[['orig.ident']]), names(sample.number.vec)),
              '\nin aggregated.csv only: ',
              setdiff(names(sample.number.vec), unique(tmp[['orig.ident']]))
              )

  )
}

# if list:
if(class(seu) == "list") {
  # 1. Rename Cells
  seu =  lapply(seu, function(seu_obj) rename.cells.seu(seu = seu_obj, rename.vec =  sample.number.vec))
} else {
  # Rename Cells
  seu =  rename.cells.seu(seu = seu, rename.vec =  sample.number.vec)
}

## 2. Save object/list 
saveRDS(seu, opt$output_seurat)

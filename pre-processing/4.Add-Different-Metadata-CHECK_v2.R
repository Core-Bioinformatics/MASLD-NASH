#!/usr/bin/env Rscript


## Add  Metadata to Seurat object or object list from the NAFLD project


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
    c('-m', "--input_metadata"),
    action="store",
    default=NA,
    type='character',
    help='Path to csv file of sample metadata.'
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

# load the seurat obj
#seu = readRDS('/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Pre-processed_data/2.Aggr_samples/2.Aggr-ALL/4.Aggr-52-samples-Aug2021/Aggr_52_samples_Aug_2021_renamed.rds')
seu = readRDS(opt$input_seurat)

# load the metadata info from table
#dt <- fread('/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/4.Analysis/1.Preprocess_data/1.Pre-process-data/Scripts/Disease status summary 63 samples updated__25_08_21_AC_cleaned.csv')
dt <- fread(opt$input_metadata)
# make the sample id, and get the columns with the metadata
meta.cols <- names(dt)[!(names(dt) %in% c('SLX', 'Sequencing.ID'))] 
dt$id <- paste0(dt$SLX, '-', dt$Sequencing.ID)
id.col <- 'id'
keep.cols <- c(id.col, meta.cols)
dt <- dt[, ..keep.cols]


add.metadata.from.table <- function(seuObj, metadata, idcol, meta.cols) {
  
  # build metadata table cell per row in same order as seuObj
  soM <- seuObj@meta.data
  soM$order <- 1:nrow(soM)
  tmp.mets <- merge(soM, metadata, all.x=T, by.x='orig.ident', by.y=idcol)
  tmp.mets <- tmp.mets[order(tmp.mets$order),]
  
  if (nrow(tmp.mets) != nrow(soM)) {
    stop('missing IDs in input_metadata table, which are in input_seurat object!')
  }
  
  # add meta.cols columns to the seuratObject
  curr.col <- meta.cols[1]
  for (curr.col in meta.cols) {
    seuObj <- AddMetaData(object=seuObj, metadata=tmp.mets[[curr.col]], col.name=curr.col)
  }
  
  return(seuObj)
}

# actually add the Annotation
# if list:
if(class(seu) == "list") {
  seu =  lapply(seu, add.metadata.from.table, metadata=dt, idcol=id.col, meta.cols=meta.cols)
} else {
  seu = add.metadata.from.table(seuObj=seu, metadata=dt, idcol=id.col, meta.cols=meta.cols)
}


## 2. Save object/list 
saveRDS(seu, opt$output_seurat)
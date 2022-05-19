
# This script runs CellPhoneDB

# Args
counts_mat=$1
meta_data=$2
#n_cells_subsample=$3
out_path=$3

echo 'counts_mat='${counts_mat}
echo 'meta_data='${meta_data}
echo 'out_path='${out_path}

# Notes:
## method can be one of: 'analysis' , 'statistical_analysis'
## counts-data  #[ensembl | gene_name | hgnc_symbol] Type of gene identifiers in the counts data

## 1. Run method with statistical analysis
# /usr/local/anaconda3/envs/cpdb/bin/cellphonedb
/usr/local/bin/cellphonedb method statistical_analysis $meta_data $counts_mat \
  --counts-data hgnc_symbol \
  --threads=8 \
  --threshold=0.1 \
  --output-path $out_path
  #--subsampling \
  #--subsampling-log false \
  #--subsampling-num-cells $n_cells_subsample \

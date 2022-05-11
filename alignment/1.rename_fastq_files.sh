#!/usr/bin/env bash

# Renames fastq files so naming convention is compatible with CellRanger

#dir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/0.data/3.Aug2021/SLX-20793/'
dir=$1

for p in ${dir}/*; do
  f=$(basename ${p})
  # will modify g and mv at end
  g=${f}

  # rename fastq suffix
  g=${g//fq/fastq}
  # remove middle bit of randomness

  g=${g//.HYKNCDRXY/}
  # replace sequencing lane nomenclature
  g=${g//.s_1./_S1_L001_}
  g=${g//.s_2./_S2_L002_}
  g=${g//.s_3./_S3_L003_}
  g=${g//.s_4./_S4_L004_}

  # replace read 1 or 2
  g=${g//r_1/R1_001}
  g=${g//r_2/R2_001}
  # replace index 1 or 2
  g=${g//i_1/I1_001}
  g=${g//i_2/I2_001}

  # replace the "." between SLX-ID, and SITT... identity. 
  g=${g/./-}
  echo ${g}
  mv ${dir}${f} ${dir}${g}
done

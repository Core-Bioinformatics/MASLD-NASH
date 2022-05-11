#!/bin/bash
set -e # exit on error
set -u # unassigned variabel is error

# slurm submission script for CellRanger count



out_dir='/servers/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/2.CR_output/'
# where to look for fastqs
fastq_dir='/servers/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/0.data/6.Feb2022/SLX-21462/'


# lists files in fastq_directory
# then cuts to keep everything before 1st "_"
# then takes unique of what remains. NB "sort" is required as uniq only strips adjacent duplicates!
for srr in `find ${fastq_dir} -type f -exec basename "{}" \; | cut -d'_' -f -1 | sort | uniq`
do
  echo 'submitting job for: '${srr}
  sbatch ./CR_count.sh ${out_dir} ${srr} ${fastq_dir} ${srr}
done

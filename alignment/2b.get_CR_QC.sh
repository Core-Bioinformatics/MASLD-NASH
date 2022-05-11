#!/usr/bin/env bash
set -eu

# copy the CR html output into a single QC directory for convenience
# renamed with name of CR output directory

# which ones to plot controlled by CR_runs csv file

# output plotted as CR_QC_summary_all.d

CR_dir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/2.CR_output/'
QC_dir='/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/3.CR_pre_analysis/CR_QC/'
CR_runs='./CR_QC_to_plot_all.csv' # all csv runs want plotted
passed_CR_runs='./CR_QC_to_plot_passed_QC.csv' # the runs which decide have passed QC

if [ ! -e ${QC_dir}html ]; then
  mkdir ${QC_dir}html
fi
if [ ! -e ${QC_dir}csv ]; then
  mkdir ${QC_dir}/csv
fi

# itrate over directories
for dir in ${CR_dir}*/; do
  outname=$(basename ${dir})
  echo ${outname}
  if [ -e ${dir}/outs/web_summary.html ]
    then cp ${dir}/outs/web_summary.html ${QC_dir}/html/${outname}.html
  fi

  if [ -e ${dir}/outs/metrics_summary.csv ]; then
    cp ${dir}/outs/metrics_summary.csv ${QC_dir}/csv/${outname}.csv
  fi
done

Rscript ./plot_CR_QC.R \
--input_csv_dir ${QC_dir}"/csv/" \
-o ${QC_dir}"CR_QC_summary.pdf" \
--input_files ${CR_runs} \
--passed_QC_csv ${passed_CR_runs}

#!/bin/bash

#SBATCH -p cccp # Partition to submit to
#SBATCH --qos=intr
#SBATCH --time=2-00:00:00 # set max time d-hh:mm:ss
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH --ntasks=1 --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=awc30@cam.ac.uk
#SBATCH -e /servers/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/0.data/slurm_out/job_%j.err
#SBATCH -o /servers/sutherland-scratch/rc845/Vallier/NAFLD-Gribben/0.data/slurm_out/job_%j.out

# args
CR_out_dir=$1 # CR output directory (2.CR_out/)
id=$2 # Output subdirectory (will be made WITHIN 2.CR_out/)
fastqs=$3 # FASTQ directory
sample_name=$4 # Fastq file prefixes within the FASTQ directory

transcriptome=/opt/bio-shares/bioinf-facility/genomes/CellRanger_references/refdata-cellranger-GRCh38-5.0.0


# CellRanger output is made in the current working directory
cd ${CR_out_dir}

echo "id: $id"
echo "fasqs: $fastqs"
echo "sample_name: $sample_name"
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

# Run CellRanger
/home/USSR/rc845/software/cellranger-5.0.0/cellranger count --id=$id \
	--fastqs=$fastqs \
	--sample=$sample_name \
	--transcriptome=$transcriptome \
	--include-introns \
 	--localmem=50 \
	--localcores=24
	##--no-bam \
	##--jobmode /servers/bigscratch/facility-scratch/rc845/Scripts/CellRanger/slurm.template

#!/bin/bash
#SBATCH --job-name=fastq_filter
#SBATCH   --output=fastq_filter.%j.out
#SBATCH    --error=fastq_filter.%j.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p dev
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > fastq_filter.$SLURM_JOBID.sh
ml load anaconda
source activate pgenlib
if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

cutoff=25

fastqs=$PI_HOME/data/nanopore-wgs-consortium/*.fq.gz
fastq_filtered=$PI_SCRATCH/projects/nanopore/ont.${cutoff}k.fq.gz
project_root=$HOME/projects/nanopore

zcat $fastqs | head -n100 |
    parallel --no-notice --jobs=8 --pipe -N4 \
	     "${project_root}/src/fastq_filter.sh -c ${cutoff} -k" | \
    gzip -9 - > ${fastq_filtered}

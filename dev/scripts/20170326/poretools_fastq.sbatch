#!/bin/bash

# This script extracts fastq file from fast5 file 

#SBATCH --job-name=poretools_fastq
#SBATCH   --output=poretools_fastq.%j.out
#SBATCH    --error=poretools_fastq.%j.err
#SBATCH --time=1-0:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > poretools_fastq.${SLURM_JOBID}.sh
ml load anaconda
source activate poretools-69ef61d
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

MEMORY=4000
THREADS=1

DATA_DIR_ROOT=$PI_SCRATCH/data/NA12878-fast5
DATA_DIR=${DATA_DIR_ROOT}'/rel3/${chr}/${fn}/'
#DATA_DIR="."

FASTQ_FILE="${DATA_DIR_ROOT}/poretools_fastq.${SLURM_JOBID}.fq.gz"

export LD_LIBRARY_PATH=/opt/rh/python27/root/usr/lib64/:$LD_LIBRARY_PATH
export HDF5_DISABLE_VERSION_CHECK=1
# the latest dev. version of poretools is built on hdf5 1.8.17

ml load hdf5/1.8.16

poretools fastq \
	  -q ${DATA_DIR} \
    | gzip - > ${FASTQ_FILE}

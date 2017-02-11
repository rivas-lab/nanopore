#!/bin/bash
#SBATCH --job-name=fast5dl
#SBATCH   --output=fast5dl.%j.out
#SBATCH    --error=fast5dl.%j.err
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > fast5dl.${SLURM_JOBID}.sh
#ml load anaconda
#source activate pgenlib
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

DATA_DIR=$PI_SCRATCH/data/NA12878-fast5/

MEMORY=4000
THREADS=1

FILE_BASE="http://s3.amazonaws.com/nanopore-human-wgs"

CHR="20"
PART="05"
FILE="${FILE_BASE}/rel3-fast5-chr${CHR}.part${PART}.tar"
wget ${FILE} -P ${DATA_DIR}

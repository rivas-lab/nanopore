#!/bin/bash
#SBATCH --job-name=fast5gzip
#SBATCH   --output=fast5gzip.%j.out
#SBATCH    --error=fast5gzip.%j.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > fast5gzip.${SLURM_JOBID}.sh
#ml load anaconda
#source activate pgenlib
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

DATA_DIR=$PI_SCRATCH/data/NA12878-fast5/

MEMORY=1000
THREADS=1

cd $DATA_DIR

CHR="20"

for PART in "05" "04" "03" "02" "01"; do
    print $PART
    gzip --best "rel3-fast5-chr${CHR}.part${PART}.tar"
done

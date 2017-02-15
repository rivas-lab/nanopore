#!/bin/bash
#SBATCH --job-name=fast5untar
#SBATCH   --output=fast5untar.%j.out
#SBATCH    --error=fast5untar.%j.err
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > fast5untar.${SLURM_JOBID}.sh
#ml load anaconda
#source activate pgenlib
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

DATA_DIR=$PI_SCRATCH/data/NA12878-fast5/

MEMORY=4000
THREADS=1

cd $DATA_DIR

CHR="20"

for PART in "01" "02" "03" "04"; do
    tar -xzf "rel3-fast5-chr${CHR}.part${PART}.tar"
done

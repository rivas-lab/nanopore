#!/bin/bash
#SBATCH --job-name=poretools_time
#SBATCH   --output=poretools_time.%j.out
#SBATCH    --error=poretools_time.%j.err
#SBATCH --time=1-0:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > poretools_time.${SLURM_JOBID}.sh
ml load anaconda
source activate poretools-69ef61d
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

DATA_DIR_ROOT=$PI_SCRATCH/data/NA12878-fast5
DATA_DIR=${DATA_DIR_ROOT}'/rel3/${chr}/${fn}/'
#DATA_DIR="."

STATS_FILE="${DATA_DIR_ROOT}/poretools_time.${SLURM_JOBID}.stats"

MEMORY=4000
THREADS=1

export LD_LIBRARY_PATH=/opt/rh/python27/root/usr/lib64/:$LD_LIBRARY_PATH
ml load hdf5/1.8.16

export HDF5_DISABLE_VERSION_CHECK=1
# the latest dev. version of poretools is built on hdf5 1.8.17

poretools times \
	  -q ${DATA_DIR} \
	  > ${STATS_FILE}

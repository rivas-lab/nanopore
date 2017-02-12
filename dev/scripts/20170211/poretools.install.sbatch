#!/bin/bash
#SBATCH --job-name=poretools.install
#SBATCH   --output=poretools.install.%j.out
#SBATCH    --error=poretools.install.%j.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p dev
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > poretools.install.${SLURM_JOBID}.sh
ml load anaconda
source activate poretools-69ef61d
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

MEMORY=4000
THREADS=1

cd $SCRATCH/poretools/69ef61d/poretools/
python2 $SCRATCH/poretools/69ef61d/poretools/setup.py install --user

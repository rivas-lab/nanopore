#!/bin/bash

# This script downloads Platinum Genome Dataset

#SBATCH --job-name=platinum-dl.dl
#SBATCH   --output=platinum-dl.dl.out
#SBATCH    --error=platinum-dl.dl.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#################
cat $0 > platinum-dl.${SLURM_JOBID}.sh

DIR=${PI_HOME}/data/PlatinumGenomes-2016v1-20170326
cd ${DIR}

wget -r -np -l 0 ftp://ussd-ftp.illumina.com/2016-1.0/

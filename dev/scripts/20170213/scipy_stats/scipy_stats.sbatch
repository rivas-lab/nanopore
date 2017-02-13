#!/bin/bash
#SBATCH --job-name=scipy_stats
#SBATCH   --output=scipy_stats.%j.out
#SBATCH    --error=scipy_stats.%j.err
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH -p dev
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > scipy_stats.${SLURM_JOBID}.sh

ml load anaconda

NAME_OF_ENV=debug-scipy_stats-${SLURM_JOBID}

conda create --yes --name ${NAME_OF_ENV} python=2.7 scipy 

source activate ${NAME_OF_ENV}

python --version

python <<EOF
from scipy import stats
EOF


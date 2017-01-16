#!/bin/bash
#SBATCH --job-name=nanopolish
#SBATCH   --output=nanopolish.out
#SBATCH    --error=nanopolish.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=8000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
target="/home/ytanigaw/data/nanopolish"

module load gcc/5.3.0
module load hdf5/1.8.16

cd ${target}
make

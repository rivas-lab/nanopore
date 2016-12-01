#!/bin/bash
#SBATCH --job-name=giab.na12878.hist
#SBATCH   --output=giab.na12878.hist.out
#SBATCH    --error=giab.na12878.hist.err
#SBATCH --time=0-10:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

DIR=${PI_HOME}/data/giab/NA12878_PacBio_MtSinai
STATS=${DIR}/sorted_final_merged.bam.stats.head
PROG=/home/ytanigaw/projects/nanopore/src/stats2hist.py

module load anaconda
source activate pgenlib

python2 ${PROG} ${STATS}

cp ${STATS}.hist.png .

git checkout plot
git add ./*.png
git commit -m "add plot of ${STATS}"
git push origin plot
git checkout master

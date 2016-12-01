#!/bin/bash
#SBATCH --job-name=giab.na12878.stats
#SBATCH   --output=giab.na12878.stats.out
#SBATCH    --error=giab.na12878.stats.err
#SBATCH --time=0-1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

DIR=${PI_HOME}/data/giab/NA12878_PacBio_MtSinai
cd ${DIR}

samtools view sorted_final_merged.bam| awk '{print $3, $5, length($10)}' > sorted_final_merged.bam.stats

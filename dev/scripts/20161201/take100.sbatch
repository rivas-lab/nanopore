#!/bin/bash
#SBATCH --job-name=sample_reads
#SBATCH   --output=sample_reads.out
#SBATCH    --error=sample_reads.err
#SBATCH --time=0-10:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

CHR="chr1"

DIR=${PI_HOME}/data/giab/NA12878_PacBio_MtSinai
BAM=${DIR}/sorted_final_merged.15kb
OUT=${BAM}.${CHR}.100


samtools view ${BAM}.bam \
    |awk -v chr=${CHR} '{if($3==chr && $9 > 15000){print $0}}' \
    |head -n100 > ${OUT}.sam

#!/bin/bash
#SBATCH --job-name=giab.na12878.extract
#SBATCH   --output=giab.na12878.extract.out
#SBATCH    --error=giab.na12878.extract.err
#SBATCH --time=0-10:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=2000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

DIR=${PI_HOME}/data/giab/NA12878_PacBio_MtSinai
BAM=${DIR}/sorted_final_merged
OUT=${BAM}.10kb

# process header lines
samtools view -H ${BAM}.bam > ${OUT}.sam

# extract reads
samtools view ${BAM}.bam \
    | grep -v ^@ \
    | awk '{if(length($10) > 10000){print $0}}' >> ${OUT}.sam

# convert to bam
samtools view -Sb ${OUT}.sam > ${OUT}.bam

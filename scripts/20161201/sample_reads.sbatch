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

DIR=${PI_HOME}/data/giab/NA12878_PacBio_MtSinai
BAM=${DIR}/sorted_final_merged
OUT=${BAM}.sample
PROG=/home/ytanigaw/projects/nanopore/src/sample_reads.pl

samtools view ${BAM}.bam \
    | ${PROG} 1804275 100 20161201 > ${OUT}.sam

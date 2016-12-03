#!/bin/bash
#SBATCH --job-name=pacbio_filter
#SBATCH   --output=pacbio_fileter.out
#SBATCH    --error=pacbio_filter.err
#SBATCH --time=2-0:0:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=63000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

ml load anaconda
source activate pgenlib

COUNT=${HOME}/projects/nanopore/src/count.py
BAM=/share/PI/mrivas/data/giab/NA12878_PacBio_MtSinai/sorted_final_merged.15kb
REF=/share/PI/mrivas/data/hg19/hg19.fa

cat <(samtools view -H ${BAM}.bam) \
    <(paste <(python2 ${COUNT} -i <(samtools view ${BAM}.bam) ${REF} \
		     | awk 'NR > 1 {print $2 / ($1 + $2)}') \
	    <(samtools view ${BAM}.bam) \
    | awk 'BEGIN{OFS="\t"} {if($1 <= 0.1) {print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}') \
    | samtools view -Sb - > ${BAM}.filter.bam

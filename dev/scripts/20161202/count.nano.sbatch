#!/bin/bash
#SBATCH --job-name=count.nano
#SBATCH   --output=count.nano.out
#SBATCH    --error=count.nano.err
#SBATCH --time=12:0:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

ml load anaconda
source activate pgenlib

COUNT=${HOME}/projects/nanopore/src/count.py
#BAM=${HOME}/data/nanopore/20161008_wgs_caucasian_48hr.10k.bwa.mapq60.20kb.chr11
#BAM=/share/PI/mrivas/data/giab/NA12878_PacBio_MtSinai/sorted_final_merged.10kb.filter
data_loc="/home/ytanigaw/data/nanopore-wgs-consortium"
BAM="${data_loc}/nanopore-wgs.25000.bam"

REF=/share/PI/mrivas/data/hg19/hg19.fa

python2 ${COUNT} -i <(samtools view ${BAM}.bam) ${REF} > ${BAM}.count

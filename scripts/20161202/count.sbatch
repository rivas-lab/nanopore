#!/bin/bash
#SBATCH --job-name=count
#SBATCH   --output=count.out
#SBATCH    --error=count.err
#SBATCH --time=0:10:00
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
#BAM=${HOME}/data/nanopore/20161008_wgs_caucasian_48hr.10k.bwa.mapq60.20kb.chr11.bam
BAM=/share/PI/mrivas/data/giab/NA12878_PacBio_MtSinai/sorted_final_merged.15kb.chr1.100.bam
REF=/share/PI/mrivas/data/hg19/hg19.fa

python2 ${COUNT} -i <(samtools view ${BAM}) ${REF}

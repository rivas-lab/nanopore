#!/bin/bash
#SBATCH --job-name=count
#SBATCH   --output=count.out
#SBATCH    --error=count.err
#SBATCH --time=0:30:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

ml load anaconda
source activate pgenlib

cat /home/ytanigaw/data/nanopore/20161008_wgs_caucasian_48hr.10k.bwa.mapq60.20kb.chr11.sam.body |python2 ${HOME}/projects/nanopore/src/count.py ~/PI_HOME/data/hg19/hg19.fa

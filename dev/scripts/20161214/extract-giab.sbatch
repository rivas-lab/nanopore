#!/bin/bash
#SBATCH --job-name=extract-giab
#SBATCH   --output=extract-giab.%j.out
#SBATCH    --error=extract-giab.%j.err
#SBATCH --time=20:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

ml load anaconda
source activate pgenlib

data_loc="${PI_HOME}/data/giab/NA12878_PacBio_MtSinai"
name="sorted_final_merged.10kb"
ref="${PI_HOME}/data/hg19/hg19.fa"

python2 ${HOME}/projects/nanopore/src/filter_by_mismatch_ratio.py \
	-i <(samtools view ${data_loc}/${name}.bam) \
	-o ${data_loc}/${name}.ext.sam \
	-e ${data_loc}/${name}.err.sam \
	-l 10000 \
	-m 0.1 \
	${ref}

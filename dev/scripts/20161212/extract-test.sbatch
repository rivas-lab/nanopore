#!/bin/bash
#SBATCH --job-name=extract-test
#SBATCH   --output=extract-test.out
#SBATCH    --error=extract-test.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

ml load anaconda
source activate pgenlib

data_loc="${HOME}/data/nanopore-wgs-consortium"
name="nanopore-wgs.25000.sorted.10k.mapq50"
ref="${PI_HOME}/data/hg19/hg19.fa"

python2 ${HOME}/projects/nanopore/src/count.py \
	--extract -q 20 \
	-i <(samtools view ${data_loc}/${name}.bam) ${ref} \
	-o ${data_loc}/${name}.ext.test.sam \
	-e ${data_loc}/${name}.err.test.sam

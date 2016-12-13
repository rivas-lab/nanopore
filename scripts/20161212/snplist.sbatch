#!/bin/bash
#SBATCH --job-name=snplist
#SBATCH   --output=snplist.out
#SBATCH    --error=snplist.err
#SBATCH --time=12:00:00
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
name="nanopore-wgs.25000.sorted.10k.mapq50.ext"
BAMFILE="${data_loc}/${name}.bam"
OUTFILE="${data_loc}/${name}.snplist"

ref="${PI_HOME}/data/hg19/hg19.fa"


python2 ../../src/dump_snps.py $ref \
	--showname \
	-i <(samtools view ${BAMFILE}) \
	-o ${OUTFILE}

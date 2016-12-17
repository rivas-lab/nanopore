#!/bin/bash
#SBATCH --job-name=dbsnp
#SBATCH   --output=dbsnp.%j.out
#SBATCH    --error=dbsnp.%j.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

cat $0 1>&2

ml load anaconda
source activate pgenlib

ml load python/2.7.5

data_loc="${HOME}/data/nanopore-wgs-consortium"
name="nanopore-wgs.25000.sorted.10k.mapq50.ext"
BAMFILE="${data_loc}/${name}.bam"

DBSNP="/share/PI/mrivas/data/dbsnp/dbsnp_all_20151104.vcf.gz"
ref="${PI_HOME}/data/hg19/hg19.fa"
WINDOW=0

OUTFILE="${data_loc}/${name}.w${WINDOW}.snps"

python2 ../../src/dump_snps.py $ref \
	-i <(samtools view ${BAMFILE}) \
	-s ${DBSNP} > ${OUTFILE}

####################################################################

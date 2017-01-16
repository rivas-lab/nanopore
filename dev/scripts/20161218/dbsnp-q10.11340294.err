#!/bin/bash
#SBATCH --job-name=dbsnp-q10
#SBATCH   --output=dbsnp-q10.%j.out
#SBATCH    --error=dbsnp-q10.%j.err
#SBATCH --time=2-0:00:00
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

data_loc="${HOME}/data/nanopore-wgs-consortium"
name="nanopore-wgs.25000.sorted.10k.mapq50.ext.sorted"
BAMFILE="${data_loc}/${name}.bam"

DBSNP="/share/PI/mrivas/data/dbsnp/dbsnp_all_20151104.vcf.gz"
ref="${PI_HOME}/data/hg19/hg19.fa"
Q_CUTOFF=10

OUTFILE="${data_loc}/${name}.q${Q_CUTOFF}.snps"

python2 ../../src/dump_snps.py $ref \
	--showname \
	-q ${Q_CUTOFF} \
	-i <(samtools view ${BAMFILE}) \
	-s ${DBSNP} > ${OUTFILE}

####################################################################

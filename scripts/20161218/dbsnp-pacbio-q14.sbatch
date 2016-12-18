#!/bin/bash
#SBATCH --job-name=dbsnp-pacbio-q14
#SBATCH   --output=dbsnp-pacbio-q14.%j.out
#SBATCH    --error=dbsnp-pacbio-q14.%j.err
#SBATCH --time=7-0:00:00
#SBATCH --qos=long
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

cat $0 1>&2

ml load anaconda
source activate pgenlib

#data_loc="${HOME}/data/nanopore-wgs-consortium"
#name="nanopore-wgs.25000.sorted.10k.mapq50.ext.sorted"
#BAMFILE="${data_loc}/${name}.bam"
data_loc="/share/PI/mrivas/data/giab/NA12878_PacBio_MtSinai"
name="sorted_final_merged.10kb.ext"
SAMFILE="${data_loc}/${name}.sam"

DBSNP="/share/PI/mrivas/data/dbsnp/dbsnp_all_20151104.vcf.gz"
ref="${PI_HOME}/data/hg19/hg19.fa"
Q_CUTOFF=14

OUTFILE="${data_loc}/${name}.q${Q_CUTOFF}.snps"

python2 ../../src/dump_snps.py $ref \
	--showname \
	-q ${Q_CUTOFF} \
	-i ${SAMFILE} \
	-s ${DBSNP} > ${OUTFILE}

####################################################################

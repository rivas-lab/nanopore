#!/bin/bash
#SBATCH --job-name=dbsnp
#SBATCH   --output=dbsnp.out
#SBATCH    --error=dbsnp.err
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
OUTFILE="${data_loc}/${name}.snps"

DBSNP="/share/PI/mrivas/data/dbsnp/dbsnp_all_20160527.vcf.gz"
ref="${PI_HOME}/data/hg19/hg19.fa"
WINDOW=2

while read QUEUE; do
    echo ">${QUEUE}"
    tabix $DBSNP $QUEUE
done < <(python2 ../../src/dump_snps.py $ref \
		 -i <(samtools view ${BAMFILE}) \
		|awk '{if(NR>1){print $0}}' \
		|sed -e 's/;/\n/g' \
		|sed -e 's/:/,/g' \
		|sed -e 's/,/ /g' \
		|awk -v window=${WINDOW} '{print $1, ":", $2-window, "-", $2+window}' \
		|sed -e 's/\s//g') > ${OUTFILE}

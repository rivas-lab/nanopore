#!/bin/bash
#SBATCH --job-name=hist
#SBATCH   --output=hist.out
#SBATCH    --error=hist.err
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
name="nanopore-wgs.25000.sorted.10k.mapq50.ext"
BAMFILE="${data_loc}/${name}.bam"
OUTFILE="${data_loc}/${name}.hist"

ref="${PI_HOME}/data/hg19/hg19.fa"

python2 ../../src/counts.py \
	-i <(samtools view ${BAMFILE}) \
	-q20 $ref \
    |awk '{if(NR > 2) {print $3 + $4}}'\
    |sort -n \
    |uniq -c > ${OUTFILE}

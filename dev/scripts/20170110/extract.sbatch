#!/bin/bash
#SBATCH --job-name=extract
#SBATCH   --output=extract.%j.out
#SBATCH    --error=extract.%j.err
#SBATCH --time=0-2:00:00
#SBATCH --qos=normal
#SBATCH -p dev
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
# This script tries to extract informative fragmentss based on UUID
# 'informative' means
#  mapped fragment length >= 10kb  &&
#  mismatch ratio <= 0.1 &&
#  has 10+ SNPs with base call Q >= 14 that are also present in dbSNP
#################

cat $0 1>&2
SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
if [ ! -e $SCRATCH_TODAY ]; then mkdir -p $SCRATCH_TODAY; fi


SRC=/share/PI/mrivas/data/nanopore-wgs-consortium/nanopore-wgs.25000.sorted.10k.mapq50.ext.sorted
TARGET=$SCRATCH_TODAY/$(basename $SRC).informative

samtools view -H ${SRC}.bam > ${TARGET}.sam
while read ID; do
    samtools view ${SRC}.bam | grep ${ID} >> ${TARGET}.sam
done < <(cat ${SRC}.q14.snps| awk '{if (NR > 1) {if($4 >= 10) {print $1}}}' | sort -u )

samtools view -Sb ${TARGET}.sam > ${TARGET}.unsorted.bam
samtools sort -o ${TARGET}.bam ${TARGET}.unsorted.bam

cp ${TARGET}.bam /share/PI/mrivas/data/nanopore-wgs-consortium/

#################

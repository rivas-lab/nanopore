#!/bin/bash
#SBATCH --job-name=extract699
#SBATCH   --output=extract699.%j.out
#SBATCH    --error=extract699.%j.err
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


BAM_SRC=/share/PI/mrivas/data/nanopore-wgs-consortium/nanopore-wgs.25000.sorted.10k.mapq50.ext.sorted.bam
SAM_TAR=/scratch/users/ytanigaw/20161224/nanopore-wgs-699.sam
BAM_TAR=/scratch/users/ytanigaw/20161224/nanopore-wgs-699.bam
IDS=/scratch/users/ytanigaw/20161224/nanopore-wgs-699.ids



samtools view -H ${BAM_SRC} > ${SAM_TAR}

while read ID; do
    samtools view ${BAM_SRC} | grep ${ID} >> ${SAM_TAR}
done < $IDS

samtools view -Sb ${SAM_TAR} > ${BAM_TAR}

cp ${BAM_TAR} /share/PI/mrivas/data/nanopore-wgs-consortium/

#################

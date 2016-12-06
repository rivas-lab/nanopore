#!/bin/bash
#SBATCH --job-name=bwa.nanopore-wgs
#SBATCH   --output=bwa.nanopore-wgs.out
#SBATCH    --error=bwa.nanopore-wgs.err
#SBATCH --time=7-0:00:00
#SBATCH --qos=long
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

ref="${PI_HOME}/data/hg19/hg19.fa"
data_loc="/home/ytanigaw/data/nanopore-wgs-consortium"
FASTQZ="${data_loc}/nanopore-wgs.25000.fq.gz"

bwa mem -x ont2d -t8 ${ref} <(zcat ${FASTQZ}) \
    | samtools view -Sb - > ${data_loc}/nanopore-wgs.25000.bam

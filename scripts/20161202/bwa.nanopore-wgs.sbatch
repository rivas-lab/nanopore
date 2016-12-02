#!/bin/bash
#SBATCH --job-name=bwa.nanopore-wgs
#SBATCH   --output=bwa.nanopore-wgs.out
#SBATCH    --error=bwa.nanopore-wgs.err
#SBATCH --time=0-1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=63000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

ref="${PI_HOME}/data/hg19/hg19.fa"
data_loc="/home/ytanigaw/data/nanopore-wgs-consortium"

bwa mem -x ont2d ${ref} <(zcat ${data_loc}/*.fq.gz) > ${data_loc}/nanopore-wgs.25000.sam

#!/bin/bash
#SBATCH --job-name=zcat.nanopore-wgs
#SBATCH   --output=zcat.nanopore-wgs.out
#SBATCH    --error=zcat.nanopore-wgs.err
#SBATCH --time=0-12:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

data_loc="/home/ytanigaw/data/nanopore-wgs-consortium"

zcat ${data_loc}/*.fq.gz \
    | gzip > ${data_loc}/nanopore-wgs.25000.fq.gz

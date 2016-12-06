#!/bin/bash
#SBATCH --job-name=sort-bam
#SBATCH   --output=sort-bam.out
#SBATCH    --error=sort-bam.err
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
data_loc="/home/ytanigaw/data/nanopore-wgs-consortium"
name="nanopore-wgs.25000"

samtools sort -@ 8 -m 4G \
	 -T ${LOCAL_SCRATCH}/${name} \
	 -o ${data_loc}/${name}.sorted.bam \
	 ${data_loc}/${name}.bam

samtools index -b ${data_loc}/${name}.sorted.bam

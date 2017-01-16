#!/bin/bash
#SBATCH --job-name=filter-bam
#SBATCH   --output=filter-bam.out
#SBATCH    --error=filter-bam.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

minQ=50
minLen=10000
minLenk=`expr $minLen / 1000`

data_loc="${HOME}/data/nanopore-wgs-consortium"
name="nanopore-wgs.25000.sorted"

cat <(samtools view -H ${data_loc}/${name}.bam) \
    <(samtools view ${data_loc}/${name}.bam \
	     | awk -v minQ="$minQ" -v minLen="$minLen" \
		   'BEGIN{OFS="\t"} {if($5 >= minQ && length($10) >= minLen) {print $0}}') \
    | samtools view -Sb - > ${data_loc}/${name}.${minLenk}k.mapq${minQ}.bam

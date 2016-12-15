#!/bin/bash
#SBATCH --job-name=dbsnp_ALL_20151104.dl
#SBATCH   --output=dbsnp_ALL_20151104.dl.%j.out
#SBATCH    --error=dbsnp_ALL_20151104.dl.%j.err
#SBATCH --time=0-3:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=2000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

# Download dbSNP from FTP server
# wget -O - \
#     ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/All_20151104.vcf.gz \
#    | zcat - \
#    | awk '(! /^#/){ print "chr"$0 }(/^#/){ print $0 }' \
#	  > /home/ytanigaw/PI_HOME/data/dbsnp/dbsnp_all_20151104.vcf

# compress
bgzip \
      /home/ytanigaw/PI_HOME/data/dbsnp/dbsnp_all_20151104.vcf \
    > /home/ytanigaw/PI_HOME/data/dbsnp/dbsnp_all_20151104.vcf.gz

# index the table file
tabix -p vcf \
      /home/ytanigaw/PI_HOME/data/dbsnp/dbsnp_all_20151104.vcf.gz

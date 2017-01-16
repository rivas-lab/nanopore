#!/bin/bash
#SBATCH --job-name=plink19.recode.vcf
#SBATCH   --output=plink19.recode.vcf.out
#SBATCH    --error=plink19.recode.vcf.err
#################
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=6:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
#
#SBATCH -p normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=4000
#################


DIR=${PI_HOME}/ukbb/download
cd ${DIR}

ml load plink/1.9_b1

for CHR in 11
do
    echo "[INFO]chrom${CHR}-plink--------------------------------"
    plink --bgen chr${CHR}impv1.bgen --sample impv1.sample --allow-no-sex --from-bp 9160000 --to-bp 9260000 --recode vcf --out chr${CHR}_region1_ont --chr ${CHR}
done


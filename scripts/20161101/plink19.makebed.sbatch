#!/bin/bash
#SBATCH --job-name=plink19.makebed
#SBATCH   --output=plink19.makebed.out
#SBATCH    --error=plink19.makebed.err
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


DIR=${PI_HOME}/data/ukbb/cal
cd ${DIR}

ml load plink/1.9_b1

for CHR in $(seq 1 22) X Y MT
do
    echo "[INFO]chrom${CHR}-plink--------------------------------"
    plink --allow-no-sex --file chrom${CHR} --make-bed --out chrom${CHR}
done


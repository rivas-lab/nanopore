#!/bin/bash
#SBATCH --job-name=ukbb.imp.chr11.md5
#SBATCH   --output=ukbb.imp.chr11.md5.out
#SBATCH    --error=ukbb.imp.chr11.md5.err
#################
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=0:30:00
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

# our "happy" chromosome :)
for CHR in 11
do
    #echo "[INFO]chrom${CHR}-gfetch-----------------------------------------"
    #./gfetch imp ${CHR} -ak1669.ukbkey
    echo "[INFO]chrom${CHR}-ukb_md5----------------------------------------"
   ./ukb_md5 chr${CHR}impv1.bgen
done

# # Reference:
# - http://biobank.ctsu.ox.ac.uk/showcase/exinfo.cgi?src=AccessingGeneticData
# - https://ibg.colorado.edu/mediawiki/index.php/UK_Biobank/Downloading_the_data

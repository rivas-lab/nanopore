#!/bin/bash
#SBATCH --job-name=na12878dl
#SBATCH   --output=na12878.out
#SBATCH    --error=na12878.err
#################
#time you think you need; default is 2 hours
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
#SBATCH --mem=2000
#################

HEAD="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf"

BED="${HEAD}.bed"
VCF="${HEAD}.vcf.gz"
TBI="${HEAD}.vcf.gz.tbi"
MD5="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NISTv3.3.md5"


DIR="${PI_HOME}/data/NA12878_GIAB_v3.3/"

mkdir -p ${DIR}

wget ${BED} -P ${DIR} 
wget ${VCF} -P ${DIR}
wget ${TBI} -P ${DIR}
wget ${MD5} -P ${DIR}

cd ${DIR}
md5sum -c NISTv3.3.md5

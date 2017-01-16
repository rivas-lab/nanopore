#!/bin/bash
#SBATCH --job-name=giab.na12878.dl
#SBATCH   --output=giab.na12878.dl.out
#SBATCH    --error=giab.na12878.dl.err
#################
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=2-0:00:00
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

DIR=${PI_HOME}/data/giab/
TARGET="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/"

wget -r -l1 -nc -nd -p ${DIR} ${TARGET}


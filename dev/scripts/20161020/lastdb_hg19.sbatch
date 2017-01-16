#!/bin/bash
#SBATCH --job-name=lastdb_hg19
#SBATCH --output=lastdb_hg19.out
#SBATCH  --error=lastdb_hg19.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=2:00:00
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
#SBATCH --mem=32000
#################

lastdb -uNEAR -R01 \
       ${PI_HOME}/data/hg19/hg19.lastdb \
       ${PI_HOME}/data/hg19/*.fa

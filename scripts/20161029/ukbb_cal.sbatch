#!/bin/bash
#SBATCH --job-name=ukbb_cal
#SBATCH   --output=ukbb_cal.out
#SBATCH    --error=ukbb_cal.err
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


DIR=${PI_HOME}/ukbb/download
cd ${DIR}

for CHR in 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26
do
    echo "${CHR}------------------------------------------"
    ./gfetch cal ${CHR} -ak1669.ukbkey
done

# # Reference:
# - http://biobank.ctsu.ox.ac.uk/showcase/exinfo.cgi?src=AccessingGeneticData
# - https://ibg.colorado.edu/mediawiki/index.php/UK_Biobank/Downloading_the_data

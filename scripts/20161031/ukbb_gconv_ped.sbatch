#!/bin/bash
#SBATCH --job-name=ukbb_gconv_ped
#SBATCH   --output=ukbb_gconv_ped.out
#SBATCH    --error=ukbb_gconv_ped.err
#################
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=4:00:00
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

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22 X Y MT
do
    echo "${CHR}------------------------------------------"
    ./gconv chrom${CHR}.cal chrom${CHR} ped
done

# # Reference:
# - http://biobank.ctsu.ox.ac.uk/showcase/exinfo.cgi?src=AccessingGeneticData
# - https://ibg.colorado.edu/mediawiki/index.php/UK_Biobank/Downloading_the_data

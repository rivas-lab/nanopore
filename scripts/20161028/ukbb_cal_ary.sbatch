#!/bin/bash
#SBATCH --job-name=ukbbcal%a
#SBATCH   --output=ukbbcal%a.out
#SBATCH    --error=ukbbcal%a.err
#################
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=12:00:00
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

####
#submit an array job
#SBATCH --array=1-24

DIR=${PI_HOME}/ukbb/download

CHR=${SLURM_ARRAY_TASK_ID}

cd ${DIR}
./gfetch cal ${CHR} -ak1669.ukbkey


# # Reference:
# - http://biobank.ctsu.ox.ac.uk/showcase/exinfo.cgi?src=AccessingGeneticData
# - https://ibg.colorado.edu/mediawiki/index.php/UK_Biobank/Downloading_the_data
# - http://slurm.schedmd.com/job_array.html

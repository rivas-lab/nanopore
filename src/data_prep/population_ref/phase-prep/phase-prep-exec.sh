#!/bin/bash

#SBATCH --job-name=phase-prep-exec
#SBATCH   --output=phase-prep-exec.%j.out
#SBATCH    --error=phase-prep-exec.%j.err
#SBATCH --time=1-0:00:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --cores=10
#SBATCH --mem=100g
#SBATCH --mail-type=END,FAIL
#################
set -beEu -o pipefail
cat $0 > phase-prep-exec.${SLURM_JOBID}.sh
#################

input="../../../../private_data/input/population_ref/chr20-geno"
output="../../../../sandbox/data_prep/population_ref/chr20-geno.bcf"

memory=90000
threads=10

script=/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/population_ref/phase-prep/phase-prep.sh

bash $script $input $output $memory $threads


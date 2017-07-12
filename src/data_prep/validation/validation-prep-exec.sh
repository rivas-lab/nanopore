#!/bin/bash

#SBATCH --job-name=validation-prep-exec
#SBATCH   --output=validation-prep-exec.%j.out
#SBATCH    --error=validation-prep-exec.%j.err
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --cores=2
#SBATCH --mem=16g
#SBATCH --mail-type=END,FAIL
#################
set -beEu -o pipefail
cat $0 > validation-prep-exec.${SLURM_JOBID}.sh
#################

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load ytanigaw-sh2

input="$nanopore/public_data/input/validation/NA12878.vcf.gz"
output="$nanopore/public_data/intermediate/validation/NA12878-snps"

memory=16000
threads=2

script=/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/validation/validation-prep.sh

bash $script $input $output $memory $threads


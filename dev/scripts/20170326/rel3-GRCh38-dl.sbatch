#!/bin/bash

# This script downloads mapped read on GRCh38 of Nanopore consortium rel3 data

#SBATCH --job-name=rel3-GRCh38-dl.dl
#SBATCH   --output=rel3-GRCh38-dl.dl.out
#SBATCH    --error=rel3-GRCh38-dl.dl.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > rel3-GRCh38-dl.${SLURM_JOBID}.sh

DIR=${PI_SCRATCH}/data/nanopore-wgs-consortium-rel3
cd ${DIR}

for CHR in $(seq 1 22) X Y M ; do
    FILE=http://s3.amazonaws.com/nanopore-human-wgs/chr${CHR}.sorted.bam
    echo $FILE >&2
    #wget $FILE
    wget $FILE.bai
done

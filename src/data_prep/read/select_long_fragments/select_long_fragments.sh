#!/bin/bash

#SBATCH --job-name=remap_from_hg38_to_hg19
#SBATCH   --output=remap_from_hg38_to_hg19.%j.out
#SBATCH    --error=remap_from_hg38_to_hg19.%j.err
#SBATCH --time=1-0:00:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --mem=2000
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-448
#################
cat $0 > remap_from_hg38_to_hg19.${SLURM_JOBID}.sh
#################

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load samtools

ary_lookup="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/read/select_long_fragments/splithg19.tsv"

###
inBam=$( cat $ary_lookup | awk -v id=$SLURM_ARRAY_TASK_ID '($1 == id){print $2}' )
outDir=$(dirname $(dirname $inBam))/$(basename $(dirname $inBam)).long
outBam=$outDir/$(basename $inBam)

echo "[task $SLURM_ARRAY_TASK_ID] input: $inBam"
echo "[task $SLURM_ARRAY_TASK_ID] out:   $outBam"

if [ ! -d $(dirname $outBam) ] ; then mkdir -p $(dirname $outBam) ; fi

bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/read/select_long_fragments/select_long_fragments-sub.sh $inBam $outBam


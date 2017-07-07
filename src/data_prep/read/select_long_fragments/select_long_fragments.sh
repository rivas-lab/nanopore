#!/bin/bash

#SBATCH --job-name=select_long_fragments
#SBATCH   --output=select_long_fragments.%j.out
#SBATCH    --error=select_long_fragments.%j.err
#SBATCH --time=0-1:00:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --mem=2000
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-448
#################
cat $0 > select_long_fragments.${SLURM_JOBID}.sh
#################

ary_lookup="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/read/select_long_fragments/splithg19.tsv"

###
inBam=$( cat $ary_lookup | awk -v id=$SLURM_ARRAY_TASK_ID '($1 == id){print $2}' )
outDir=$(dirname $(dirname $inBam))/$(basename $(dirname $inBam)).long2
outBam=$outDir/$(basename $inBam)

echo "[task $SLURM_ARRAY_TASK_ID] input: $inBam"
echo "[task $SLURM_ARRAY_TASK_ID] out:   $outBam"

if [ ! -d $(dirname $outBam) ] ; then mkdir -p $(dirname $outBam) ; fi

bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/read/select_long_fragments/select_long_fragments-sub.sh $inBam $outBam


#!/bin/bash

#SBATCH --job-name=remap_from_hg38_to_hg19
#SBATCH   --output=remap_from_hg38_to_hg19.%j.out
#SBATCH    --error=remap_from_hg38_to_hg19.%j.err
#SBATCH --time=1-0:00:00
#SBATCH --qos=normal
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --cores=10
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-448
#################
cat $0 > remap_from_hg38_to_hg19.${SLURM_JOBID}.sh
#################

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load samtools

n_threads=10
mem=64000

ref="/oak/stanford/groups/mrivas/public_data/hg19/chr20.fa"
ary_lookup="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/read/remap_from_hg38_to_hg19/splitSam.tsv"

###
inSam=$( cat $ary_lookup | awk -v id=$SLURM_ARRAY_TASK_ID '($1 == id){print $2}' )
outDir=$(dirname $(dirname $inSam))/$(basename $(dirname $inSam)).hg19
outBam=$outDir/$(basename ${inSam%.sam}).bam

tmpDir=$(mktemp -d); echo "[$0] tmpDir is $tmpDir"
handler_exit () { rm -rf $tmpDir ; }
trap handler_exit EXIT


tmpBam1=$(mktemp -p $tmpDir)
tmpBam2=$(mktemp -p $tmpDir)

echo "[task $SLURM_ARRAY_TASK_ID] input: $inSam"
echo "[task $SLURM_ARRAY_TASK_ID] tmp1:  $tmpBam1"
echo "[task $SLURM_ARRAY_TASK_ID] tmp2:  $tmpBam2"
echo "[task $SLURM_ARRAY_TASK_ID] out:   $outBam"

echo "[$0] mapping..."
samtools fastq $inSam  \
        | bwa mem -x ont2d -t $n_threads $ref - \
        | samtools view -Sb - > $tmpBam1

echo "[$0] sorting..."
samtools sort -l 9 -@ $n_threads -m ${mem}M -o $tmpBam2 $tmpBam1

echo "[$0] moving results file"
if [ ! -d $(dirname $outBam) ] ; then mkdir -p $(dirname $outBam) ; fi
mv $tmpBam2 $outBam

#################

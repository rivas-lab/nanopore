#!/bin/bash
set -beEu -o pipefail


# input 

input=$1
output=$2

memory=$3
threads=$4

####

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load plink2

which plink2
plink2 --version

####

echo "[$0] running PLINK2 to export to a pgen/bim/fam file"
plink2 --vcf $input \
	--memory $memory --threads $threads \
	--snps-only just-acgt \
	--chr 20 \
	--out $output \
	--make-bpgen


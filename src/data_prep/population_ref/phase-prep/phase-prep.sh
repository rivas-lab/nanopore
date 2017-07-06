#!/bin/bash
set -beEu -o pipefail


# input 

input=$1
output=$2

memory=$3
threads=$4

####

tmpDir=$(mktemp -d ); echo "[$0] tmpDir is $tmpDir"
handler_exit () { rm -rf $tmpDir ; }
trap handler_exit EXIT

####

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load plink2
module load bcftools

which plink2
plink2 --version
which bcftools
bcftools --version

####

tmp_vcf=$tmpDir/tmp-vcf.vcf
tmp_bcf=$tmpDir/tmp-bcf.bcf

plink2 --bfile $input --export vcf --memory $memory --threads $threads --out ${tmp_vcf%.vcf}

bcftools convert --output-type b --threads $threads --output $tmp_bcf $tmp_vcf

mv $tmp_bcf $output


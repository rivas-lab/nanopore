#!/bin/bash
set -beEu -o pipefail


# input 

input=$1
output=$2
keep=$3

memory=$4
threads=$5

####

tmpDir=$(mktemp -d ); echo "[$0] tmpDir is $tmpDir"
handler_exit () { rm -rf $tmpDir ; }
trap handler_exit EXIT

####

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load plink
module load bcftools

which plink
plink --version
which bcftools
bcftools --version

####

tmp_vcf=$tmpDir/tmp-vcf.vcf
tmp_bcf=$tmpDir/tmp-bcf.bcf


echo "[$0] running PLINK2 to export to a vcf file"
plink --bfile $input \
	--memory $memory --threads $threads \
	--keep $keep \
	--geno --hwe 1e-10 midp --maf 0.005 --biallelic-only strict \
	--snps-only just-acgt \
	--out ${tmp_vcf%.vcf} --export vcf

echo "[$0] running bcftools to convert to a bcf file"
bcftools convert --output-type b --threads $threads --output $tmp_bcf $tmp_vcf

mv $tmp_bcf $output


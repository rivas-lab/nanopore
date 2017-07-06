#!/bin/bash
set -beEu -o pipefail

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load samtools

bam=$1
outFile=$2
minLen=10000

outTmp=$(mktemp)
#echo $bam $outTmp $outFile

samtools view -h $bam \
	| awk -v minLen=$minLen '((substr($0, 1, 1) == "@") || length($10) >= minLen){print $0}' \
	| samtools view -bS - > $outTmp
mv $outTmp $outFile


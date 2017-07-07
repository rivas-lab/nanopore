#!/bin/bash
set -beEu -o pipefail

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
module load samtools

bam=$1
outFile=$2
minLen=10000
minMapQ=50

outTmp=$(mktemp)

echo "[$0] $bam" >&2
echo "[$0] $outTmp" >&2
echo "[$0] $outFile" >&2

echo "[$0] samtools view"
samtools view -q $minMapQ -h $bam \
	| awk -v minLen=$minLen '((substr($0, 1, 1) == "@") || length($10) >= minLen){print $0}' \
	| samtools view -bS - > $outTmp

echo "[$0] writing to the results file" > &2
mv $outTmp $outFile


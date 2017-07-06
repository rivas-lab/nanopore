#!/bin/bash
set -beEu -o pipefail

if [ $# -lt 2 ] ; then echo "usage: $0 <in.bam> <outDir>" >&2 ; exit 1 ; fi

inBam=$(readlink -m $1)
outDir=$2

if [ -d $outDir ] ; then mv $outDir $outDir-$(date +%Y%m%d-%H%M%S) ; fi

# make a temp dir
tmpDir=$(mktemp -d); echo "[$0] tmpDir is $tmpDir"
handler_exit () { rm -rf $tmpDir ; }
trap handler_exit EXIT


mkdir -p $tmpDir/txt $tmpDir/sam

echo "[$0] samtools version check"
which samtools
samtools --version 


# dump header
samtools view -H $inBam > $tmpDir/sam.head


echo "[$0] split the input file"
split -d -l 1000 <( samtools view $inBam ) $tmpDir/txt/$(basename $inBam)- 


find $tmpDir/txt -type f > $tmpDir/fileList
echo "[$0] there are $(cat $tmpDir/fileList | wc -l) files"


while read subTxt ; do
	baseName=$(basename $subTxt)
	cat $tmpDir/sam.head $subTxt > $tmpDir/sam/$baseName.sam
done < $tmpDir/fileList


echo "[$0] writing results to $outDir"
cp -r $tmpDir/sam/ $outDir

#!/bin/bash

fasta=$1
bim=$2
bam=$3

hap=${bam%.bam}.hap
err=${bam%.bam}.err

echo "[$0] $bam" >&2

python2 bamToHap.py --fasta $fasta --bim $bim $bam >$hap 2>$err


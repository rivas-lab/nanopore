#!/bin/bash

fasta=./data/chr
bim=./data/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes-pgen.bim
bam=./data/rel3.chr20.12500.10k.head.bam

hap=${bam%.bam}.hap
err=${bam%.bam}.err

python2 bamToHap.py --fasta $fasta --bim $bim $bam >$hap 2>$err


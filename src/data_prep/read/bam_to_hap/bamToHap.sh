#!/bin/bash

fasta=/oak/stanford/groups/mrivas/public_data/hg19/chr
bim=../../../../private_data/intermediate/population_ref/chr20-alleles.bim 
bam=../../../../sandbox/data_prep/read/remap_from_hg38_to_hg19/chr20.sorted.split.hg19.long2/chr20.sorted.bam-00.bam
hap=${bam%.bam}.hap
err=${bam%.bam}.err

python2 bamToHap.py --fasta $fasta --bim $bim $bam 2>$err | tee $hap


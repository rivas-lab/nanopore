#!/bin/sh

ml load ytanigaw-sh2
find $nanopore/sandbox/data_prep/read/remap_from_hg38_to_hg19/chr20.sorted.split.hg19.long -name "*.bam"|sort -g |awk -v OFS='\t'  '{print NR, $1}'

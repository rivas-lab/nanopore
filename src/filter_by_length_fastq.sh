#!/bin/sh

cutoff=$1

awk -v cutoff="$cutoff" 'BEGIN {OFS = "\n"} {head = $0; getline seq; getline qhead; getline qseq; if (length(seq) >= cutoff) {print head, seq, qhead, qseq}}'

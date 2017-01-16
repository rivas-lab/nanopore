#!/bin/bash
# @(#) This script takes gzipped fasta files (or stdin) and extract reads based on the length of the reads.

#-------------------
# default parameters
#-------------------
# default cutoff length is 25k
cutoff=25
base_letter="k"
# verbose mode is false
verbose=-1

#-------------------
# read cmd_args
#-------------------
while getopts ":c:kv" opt; do
    case $opt in
	c)
	    cutoff=$OPTARG
	    nshift=$(echo "${nshift:-0} + 2" | bc)
	    ;;
	k)
	    base_letter="k"
	    nshift=$(echo "${nshift:-0} + 1" | bc)
	    ;;
	v)
	    verbose=1
	    nshift=$(echo "${nshift:-0} + 1" | bc)
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    ;;
    esac
done

#-------------------
# prepare params
#-------------------
base=1
case $base_letter in
    k)
	base=1000
        ;;
    m)
	base=1000000
	;;
esac
shift $nshift
cutoff_val=$(echo "${cutoff} * ${base}" | bc)
suffix="${cutoff}${base_letter}"


#-------------------
# verbose mode
#-------------------
if [ "${verbose}" -eq 1 ]; then
    echo $(basename $0) >&2
    echo "cutoff is: $cutoff_val $suffix" >&2
    if [ -p /dev/stdin ]; then
	echo "read from stdin" >&2
    else
	echo "read from the following files: $@" >&2
    fi
fi


#-------------------
# filter fasta file with awk
#-------------------
if [ -p /dev/stdin ]; then
    cat -
else
    zcat $@
fi | awk -v cutoff="$cutoff_val" 'BEGIN {OFS = "\n"} {head = $0; getline seq; getline qhead; getline qseq; if (length(seq) >= cutoff) {print head, seq, qhead, qseq}}'

#!/bin/bash
# @(#) This script takes sam files (or stdin) and extract reads based on the length of the reads.

#-------------------
# default parameters
#-------------------
# default cutoff length is 10k
cutoff=10
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
	    nshift=$(echo "${nshift:-0} + 2" | /share/PI/mrivas/bin/bc)
	    ;;
	k)
	    base_letter="k"
	    nshift=$(echo "${nshift:-0} + 1" | /share/PI/mrivas/bin/bc)
	    ;;
	v)
	    verbose=1
	    nshift=$(echo "${nshift:-0} + 1" | /share/PI/mrivas/bin/bc)
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
cutoff_val=$(echo "${cutoff} * ${base}" | /share/PI/mrivas/bin/bc)
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
# filter bam file with awk
#-------------------
if [ -p /dev/stdin ]; then
    cat -
else
    zcat $@
fi | awk -v cutoff=$cutoff_val 'BEGIN{OFS="\t"} {if ((length($10) > cutoff) && (($2 / 2048) < 1)) {print $0}}'

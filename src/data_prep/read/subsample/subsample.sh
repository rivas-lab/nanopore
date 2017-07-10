#!/bin/sh
set -beEu -o pipefail

inHap=$1
coverage=$2
bim=$3
outHap=$4

seed=100

inHap_total_obs=$( cat $inHap | awk 'BEGIN{sum=0} {sum = sum + $8 - $7} END{print sum}' )
bim_wc=$( cat $bim | wc -l )

sampling_rate=$( python -c "print( 1.0 *  $bim_wc * $coverage / $inHap_total_obs )" )

cat $inHap | awk -v seed=$seed -v sampling_rate=$sampling_rate 'BEGIN{srand(seed)} (rand() <= sampling_rate){print $0}' > $outHap


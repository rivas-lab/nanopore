#!/bin/bash
set -beEu -o pipefail

inBam=$1

outDir=$(dirname $(dirname $inBam))/$(basename $(dirname $inBam)).long.parallel
outBam=$outDir/$(basename $inBam)

if [ ! -d $outDir ] ; then mkdir -p $outDir ; fi

bash $(dirname $(readlink -m $0))/$(basename ${0%-parallel.sh}.sh) $inBam $outBam


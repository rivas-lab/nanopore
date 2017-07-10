#!/bin/sh
set -beEu -o pipefail

export MODULEPATH=$HOME/.modules:$PI_HOME/.modules:$MODULEPATH
ml load ytanigaw-sh2

inHap=$nanopore/private_data/intermediate/read/chr20.sorted.hap
bim=$nanopore/private_data/intermediate/population_ref/chr20-alleles.bim

for coverage in 1 2 3 4 5 8 10 ; do
	outHap=$nanopore/private_data/intermediate/read/chr20.sorted.${coverage}x.hap
	echo "$outHap" >&2
	bash $(dirname $(readlink -m $0))/$(basename ${0%-exec.sh}).sh $inHap $coverage $bim $outHap
done


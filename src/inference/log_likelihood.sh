#!/bin/bash
set -beEu -o pipefail

. $(dirname $(readlink -m $0))/params.sh

python $(basename ${0%.sh}).py --hap $hap_f --bins $bins_f --key $hapkey_d  --ll $log_likelihood_d 2>&1 | tee $(basename ${0%.sh}).out


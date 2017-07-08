#!/bin/bash
set -beEu -o pipefail

. $(dirname $(readlink -m $0))/params.sh

python $(basename ${0%.sh}).py --bins $bins_f --cnt $hapcnt_d  --ll $log_likelihood_d --lp $log_posterior_d 2>&1 | tee $(basename ${0%.sh}).out


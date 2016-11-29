ml load anaconda/anaconda3
source activate pgenlib

which python
python --version
which R
R --version

cd /home/ytanigaw/syzygy-1.2.7

ENV_LOC="/home/ytanigaw/.conda/envs/pgenlib"

RLIB=${ENV_LOC}/lib/R/library/
RBIN=${ENV_LOC}/bin/R
RSCRIPT=${ENV_LOC}/bin/Rscript
PYTHON=${ENV_LOC}/bin/python
PYTHONPATH=${ENV_LOC}/lib/python2.7/site-packages/

export RLIB
export RBIN
export RSCRIPT
export PYTHON
export PYTHONPATH

make clean

./configure --prefix=$HOME --libdir=${ENV_LOC}/lib

make install

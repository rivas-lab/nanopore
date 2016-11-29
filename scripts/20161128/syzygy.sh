ml load anaconda/anaconda3
source activate pgenlib

which python
python --version
which R
R --version

cd /home/ytanigaw/syzygy-1.2.7

ENV_LOC="/home/ytanigaw/.conda/envs/pgenlib"

export RLIB=${ENV_LOC}/lib/R/library/
export RBIN=${ENV_LOC}/bin/R
export RSCRIPT=${ENV_LOC}/bin/Rscript
export PYTHON=${ENV_LOC}/bin/python
export PYTHONPATH=${ENV_LOC}/lib/python2.7/site-packages/

make clean

./configure --prefix=$HOME --libdir=${ENV_LOC}/lib

make install

ml load anaconda/anaconda3
source activate pgenlib

which python
python --version
which R
R --version

cd /home/ytanigaw/syzygy-1.2.7

RLIB=/home/ytanigaw/.conda/envs/pgenlib/lib/R/library/
RBIN=/home/ytanigaw/.conda/envs/pgenlib/bin/R
PYTHON=/home/ytanigaw/.conda/envs/pgenlib/bin/python
PYTHONPATH=/home/ytanigaw/.conda/envs/pgenlib/lib/python2.7/site-packages/

make clean

./configure --prefix=$HOME --libdir=${HOME}/.conda/envs/pgenlib/lib

make install

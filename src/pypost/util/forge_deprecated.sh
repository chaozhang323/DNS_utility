#!/bin/bash
module load openmpi/2.1.0/intel/2015.2.164
export PATH=/duan/data/yl8bc/python/Python-2.7.11:$PATH
export PYTHONPATH=/duan/data/yl8bc/python/Python-2.7.11
#export PATH=/duan/data/yl8bc/.local/bin:$PATH
module load tecplot/2017


## tecplot
tec360_lib=/share/apps/tecplot360/2017/360ex_2017r2/bin/:/share/apps/tecplot360/2017/360ex_2017r2/bin/sys/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$tec360_lib

curdir=$(pwd)
cd $( dirname $BASH_SOURCE )
python setup.py install --user
cd $curdir

echo "Hello $USER!"
echo "The script is for initialization of envrionment variables needed for Python2.7. If you see no error reported so far, you are good to proceed."

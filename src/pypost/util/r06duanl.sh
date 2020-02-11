#!/bin/bash
tec360_lib=/usr/local/tecplot/360ex_2017r3/bin/:/usr/local/tecplot/360ex_2017r3/bin/sys/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$tec360_lib

curdir=$(pwd)
cd $( dirname $BASH_SOURCE )
python3 get-pip.py --user
python3 -m pip install --user -r requirements.txt

cd util
. ../make_util
cd ..

python3 setup.py install --user
cd $curdir

echo "Hello $USER!"
echo "The script is for initialization of envrionment variables needed for Python3.5. If you see no error reported so far, you are good to proceed."


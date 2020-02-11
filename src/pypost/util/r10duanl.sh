#!/bin/bash
tec360_lib=/opt/tecplot/360ex_2017r2/bin/:/opt/tecplot/360ex_2017r2/bin/sys/

#if [[ "$LD_LIBRARY_PATH" != *"$itec360_lib"* ]]; then
#    echo "Adding tecplot360 path..."
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$tec360_lib
#fi

curdir=$(pwd)
cd $( dirname $BASH_SOURCE )
python3 setup.py install --user
cd $curdir

echo "Hello $USER!"
echo "The script is for initialization of envrionment variables needed for Python2.7. If you see no error reported so far, you are good to proceed."


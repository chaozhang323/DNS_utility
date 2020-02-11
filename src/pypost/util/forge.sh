#!/bin/bash
#!/bin/bash
#module load openmpi/2.1.1/intel/2015.2.164
module purge
module load python/3.6.4
module load tecplot/2017

## tecplot
tec360_lib=/share/apps/tecplot360/2017/360ex_2017r2/bin/:/share/apps/tecplot360/2017/360ex_2017r2/bin/sys/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$tec360_lib

curdir=$(pwd)
cd $( dirname $BASH_SOURCE )
python3 get-pip.py --user
python3 -m pip install --user --upgrade pip==9.0.3 # Broken pip10.0
pip3.6 install --user -r requirements.txt

cd util
. ../make_util
cd ..

python3 setup.py install --user
cd $curdir

echo "Hello $USER!"
echo "The script is for initialization of envrionment variables needed for Python3.5. If you see no error reported so far, you are good to proceed."

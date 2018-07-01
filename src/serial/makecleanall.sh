#!/bin/bash
echo "Usage: makecleanall [delete data files(yes/no(default))]"
workdir=`pwd`
deletedata=no
if [ $# -ge 1 ]
  then
  deletedata=$1
  if [ $deletedata == "yes" ]
    then
    echo "Data files will be delete in all sub directories"
    echo "Do you want to continue? yes/no"
    read ACTION
    if [ $ACTION != "yes" ]
      then
      exit 0
    fi
  fi
fi
for arg in $(find -type d)
do
  echo "Entering directory "$arg
  cd $arg
  make clean
  if [ $deletedata == "yes" ]
    then
    echo "deleting all .dat files"
    rm -f *.dat
  fi
  cd $workdir
done

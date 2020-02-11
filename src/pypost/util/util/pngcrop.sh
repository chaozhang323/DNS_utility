#!/bin/bash

if [ $# -le 0 ]; then
    echo
    echo "Usage: $(basename $0) file1.png [file2.png ...]"
    echo
    echo "  This script uses gimp to autocrop PNG files and"
    echo "  save them to PNG format.  You must have"
    echo "  crop-png.scm installed in your gimp "
    echo "  scripts directory."
    echo
    exit 1
fi

# set the filelist
files=$*

# # set the base command
# CMD="gimp -i -b "

# loop and add each file
for i in ${files[*]} ; do
  # #echo $i
  # ARGS="\"(crop-png \\\"$i\\\")\""
  # CMD="$CMD $ARGS"

  gimp -i -b "(crop-png \"$i\")" -b "(gimp-quit 0)"
done

# # add the end to quit
# TAIL="-b \"(gimp-quit 0)\""
# CMD="$CMD $TAIL"
# 
# #echo $CMD
# eval $CMD

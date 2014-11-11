#!/bin/bash

# This script is used as a forward step in optimization using Dinver.
# Main.py is name of python program which does the 
# calculation of misfit function for every given parameters.
# parameters are output of Dinver which are generated during the optimization 
# process.
# outfile1 and outfile2 are parameters and misfit data for every run.
# last modify: 11 Nov 2014 (adding explanations)

outfile1=$1
outfile2=$2



if [ ! -f $outfile1 ];
then
  touch $outfile1
fi
if [ ! -f $outfile2 ];
then
  touch $outfile2
fi


tr '\n' ' ' < parameters 
echo $(<parameters) >> $outfile1

python  Main.py > filename
 
tr '\n' ' ' < filename
echo $(<filename) >> $outfile2
echo $(<filename) > misfit


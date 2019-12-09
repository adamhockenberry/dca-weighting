#!/usr/bin/bash

#This little script just removes the last line from all the .mat outfiles in the indicated directory
#and then shunts them to a file: ./runs.log

#I subsequently rename/move this file over to the Results folder as CCMPredPy_MetaData.log

#This exists since the CCMPredPy output is not just a matrix of the couplings but instead includes the meta
#data line at the bottom of the file which made life a bit difficult without moving things around

touch runs.log
for couplings_file in ../Results/couplings_revision/*.mat; do
    tail -q -n1 $couplings_file >> runs.log
done
for couplings_file in ../Results/couplings_revision/*.mat; do
    sed -i '' -e '$ d' $couplings_file
done

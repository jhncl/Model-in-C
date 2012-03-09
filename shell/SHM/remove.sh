#!/bin/bash
list=`ls -d */`
for file in $list
do
cd "/home/b0919573/QFADatasets/SHM/""$file"
rm priors.txt
done
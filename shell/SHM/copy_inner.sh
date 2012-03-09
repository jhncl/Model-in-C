#!/bin/bash            
IFS=$'\n'
list=`ls -d */` 
for folder in $list
do
cd $folder
list_inner=`ls -d */`
for folder_inner in $list_inner
do
cp ~/QFADatasets/SHM/priors.R ~/QFADatasets/SHM/$folder$folder_inner.
done
cd ..
done
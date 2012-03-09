#!/bin/bash            
IFS=$'\n'
list=`ls -d */` 
for folder in $list
do
cd $folder
list_inner=`ls -d */`
for folder_inner in $list_inner
do
cd $folder_inner
rm para_*
cd ..
done
cd ..
done
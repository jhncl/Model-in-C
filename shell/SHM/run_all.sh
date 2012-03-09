#!/bin/bash            
IFS=$'\n'
list=`ls -d */` 
for folder in $list
do
cd $folder
inner_list=`ls -d */`
for inner_folder in $inner_list
do
cd "$inner_folder" 
make
sleep 1
./shellC.sh
cd ..
done
cd ..
done


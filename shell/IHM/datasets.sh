#!/bin/bash            
IFS=$'\n'
rm datasets.txt
cd ~/QFADatasets/SHM
list=`ls -d */` 
for folder in $list
do
echo $folder
echo $folder >> ~/QFADatasets/IHM/datasets.txt
cd $folder
inner_list=`ls -d */`
echo $inner_list
echo $inner_list >> ~/QFADatasets/IHM/datasets.txt
echo " "
echo " " >> ~/QFADatasets/IHM/datasets.txt
cd ..
done


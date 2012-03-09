#!/bin/bash            
IFS=$'\n'
rm datasets.txt
cd ~/QFADatasets/SHM
list=`ls -d */` 
for folder in $list
do
echo $folder
echo $folder >> ~/QFADatasets/JHM/datasets.txt
cd $folder
inner_list=`ls -d */`
echo $inner_list
echo $inner_list >> ~/QFADatasets/JHM/datasets.txt
echo " "
echo " " >> ~/QFADatasets/JHM/datasets.txt
cd ..
done


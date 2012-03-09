#!/bin/bash            
IFS=$'\n'
list=`ls -d */` 
for folder in $list
do
cd $folder
inner_list=`ls -d */`
for inner_folder in $inner_list
do
cp ~/2/Model-in-C/CCode/SHM/*.c ~/QFADatasets/SHM/$folder"$inner_folder".
cp ~/2/Model-in-C/CCode/SHM/shellC.sh ~/QFADatasets/SHM/$folder"$inner_folder".
cp ~/2/Model-in-C/CCode/SHM/*.h ~/QFADatasets/SHM/$folder"$inner_folder".
cp ~/2/Model-in-C/CCode/SHM/main ~/QFADatasets/SHM/$folder"$inner_folder".
cp ~/2/Model-in-C/CCode/SHM/Makefile ~/QFADatasets/SHM/$folder"$inner_folder".
cp ~/QFADatasets/SHM/priors.txt ~/QFADatasets/SHM/$folder"$inner_folder".
cd $inner_folder
make
cd ..
done
cd ..
done


#!/bin/bash            
IFS=$'\n'
folder=`echo "$1_$3"`
mkdir $folder
folder_inner=`echo "$2_$4"`
mkdir $folder/$folder_inner
cp ~/QFADatasets/SHM/$1/data.txt ~/QFADatasets/JHM/$folder/$folder_inner/dataA1.txt
cp ~/QFADatasets/SHM/$3/data.txt ~/QFADatasets/JHM/$folder/$folder_inner/dataB1.txt
cp ~/QFADatasets/SHM/priors_JHM.txt ~/QFADatasets/JHM/$folder/$folder_inner/priors.txt
cp ~/QFADatasets/JHM/M_JHM_FULL.R ~/QFADatasets/JHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/JHM/*.c ~/QFADatasets/JHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/JHM/*.h ~/QFADatasets/JHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/JHM/main ~/QFADatasets/JHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/JHM/Makefile ~/QFADatasets/JHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/JHM/shellC.sh ~/QFADatasets/JHM/$folder/$folder_inner/.
cd $folder/$folder_inner
make






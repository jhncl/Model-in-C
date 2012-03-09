#!/bin/bash            
list_txt=`ls *.txt` 
for file in $list_txt
do
name=`basename $file`
name=`echo $name | sed '$s/....$//'`
mkdir $name
mv $file $name/data.txt
cp M_SHM_FULL.R $name/M_SHM_FULL.R
done

#!/bin/bash            
list=`ls -d */` 
for file in $list
do
cp M_SHM_FULL.R "$file""M_SHM_FULL.R"
cp priors.R "$file""priors.R"
done
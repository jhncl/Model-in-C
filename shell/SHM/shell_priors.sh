#!/bin/bash
IFS=$'\n'
list=`ls -d */`
for folder in $list
do
cd "/home/b0919573/QFADatasets/SHM/""$folder"
list_inner=`ls -d */`
for folder_inner in $list_inner
do
cd $folder_inner
treat=`echo $folder_inner | sed '$s/.$//'`
name=`echo $folder_$treat`  
echo "setwd(\"/home/b0919573/QFADatasets/SHM/$folder$folder_inner\")
getwd()
library(gtools,lib=\"/home/b0919573/R/x86_64-pc-linux-gnu-library/2.10\")
library(gdata,lib=\"/home/b0919573/R/x86_64-pc-linux-gnu-library/2.10\")
library(gplots,lib=\"/home/b0919573/R/x86_64-pc-linux-gnu-library/2.10\") 
library(qfa,lib=\"/home/b0919573/R\")
library(qfaBayes,lib=\"/home/b0919573/R\")
options(warn=-1)
load(\"M_SHM_FULL_$treat.RData\")
source(\"priors.R\")" > I$name.R

echo "#!/bin/sh
hostname
date
R CMD BATCH I$name.R $name.log
die" > T$name.sh

echo "#####################

# My First Submission

#####################

Executable = T$name.sh

#Arguments = I bring you tidings of joy

Universe = vanilla

#Requirements = OpSys == “LINUX” && Arch ==“X86_64”

Error = hello_condor.err

Output = hello_condor.out

Log = hello_condor.log

transfer_input_files = I$name.R

transfer_files = ALWAYS

Queue
" > something.classad

chmod 755 T$name.sh
condor_submit something.classad
sleep 1
cd ..
done
done
#!/bin/bash
list=`ls -d */`
for folder in $list
do
cd "/home/b0919573/QFADatasets/SHM/""$folder"
name=`echo $folder | sed '$s/.$//'`
echo "
data_dir=\"$folder\"
setwd(\"/home/b0919573/QFADatasets/SHM/$folder\")
getwd()
library(qfa,lib=\"/home/b0919573/R\")
library(qfaBayes,lib=\"/home/b0919573/R\")
source(\"M_SHM_FULL.R\")" > I$name.R

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

done
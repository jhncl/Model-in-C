#!/bin/bash
list=`ls -d */`
for folder in $list
do
cd /home/b0919573/QFADatasets/JHM/$folder
list_inner=`ls -d */`
for folder_inner in $list_inner
do
cd $folder_inner
step1=`echo $folder | sed '$s/.$//'`
step2=`echo $folder_inner | sed '$s/.$//'`
name=$step1"_"$step2
TREAT=$folder_inner
TREAT=`echo $TREAT | sed '$s/.$//'`
TreatA=${TREAT/_*}
TreatB=${TREAT/*_}

echo "
TreatA=\"$TreatA\"
TreatB=\"$TreatB\"  
setwd(\"/home/b0919573/QFADatasets/JHM/$folder$folder_inner\")
getwd()
library(qfa,lib=\"/home/b0919573/R\")
library(qfaBayes,lib=\"/home/b0919573/R\")
source(\"M_JHM_FULL.R\")" > I$name.R

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
cd ..
done
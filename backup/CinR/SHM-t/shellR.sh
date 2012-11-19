#!/bin/bash

name=`echo $1 | sed '$s/..$//'`
echo "$name"
echo "#!/bin/sh
hostname
date
R CMD BATCH $1 R$name.log
date" > R$name.sh

echo "#####################

# My First Submission

#####################

Executable = R$name.sh

#Arguments = I bring you tidings of joy

Universe = vanilla

#Requirements = OpSys == “LINUX” && Arch ==“X86_64”

Error = hello_condor.err

Output = hello_condor.out

Log = hello_condor.log

transfer_input_files = $1, priors.txt, M_SHM_FULL_27.RData

transfer_files = ALWAYS

Queue
" > something.classad

chmod 755 R$name.sh
condor_submit something.classad
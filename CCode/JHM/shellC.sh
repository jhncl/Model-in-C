#!/bin/bash

echo "#!/bin/sh
hostname
rm out_C.R
make
sleep 10
date
chmod 755 main
./main>out_X_XXLC.R
date
die" > main.sh

echo "#####################

# My First Submission

#####################

Executable = main.sh

#Arguments = I bring you tidings of joy

Universe = vanilla

#Requirements = OpSys == “LINUX” && Arch ==“X86_64”

Error = hello_condor.err

Output = hello_condor.out

Log = hello_condor.log

transfer_input_files = main, Makefile, ydata.txt, xdata.txt, NoTIMEdata.txt, NoORFdata.txt, main.c, data2.c, functions.c, print.c, struct.c, main.c

transfer_files = ALWAYS

Queue
" > something.classad

chmod 755 main.sh
condor_submit something.classad
sleep 1

#!/bin/bash                                                                                                                                                  
step1="`pwd`"
data="`basename $step1`"
step2="`dirname $step1`"
step3="`basename $step2`"
data=$step3_$data
model="JHM"
burn="100000"
iter="1000"
thin="100"

name="CC_${model}_${data}_${burn}_${iter}_${thin}"


echo "#!/bin/sh
                                                                                                                                              
hostname                                                                                                                                                     
make                                                                                                                                                         
sleep 10                                                                                                                                                     
date                                                                                                                                                         
chmod 755 main                                                                                                                                               
./main $burn $iter $thin >${name}.R                                                                                                                             
date                                                                                                                                                         
die" > ${name}.sh

echo "#####################                                                                                                                                  
                                                                                                                                                             
# My First Submission                                                                                                                                        
                                                                                                                                                             
#####################                                                                                                                                        
                                                                                                                                                             
Executable = ${name}.sh                                                                                                                              
                                                                                                                                                             
#Arguments = I bring you tidings of joy                                                                                                                      
                                                                                                                                                             
Universe = vanilla                                                                                                                                           
                                                                                                                                                             
#Requirements = OpSys == “LINUX” && Arch ==“X86_64”                                                                                                          
                                                                                                                                                             
Error = hello_condor.err                                                                                                                                     
                                                                                                                                                             
Output = hello_condor.out                                                                                                                                    
                                                                                                                                                             
Log = hello_condor.log                                                                                                                                       
                                                                                                                                                             
transfer_input_files = main, Makefile, ydataA1.txt, ydataB1.txt, xdataA1.txt, xdataB1.txt,NoTIMEdataA1.txt,NoTIMEdataB1.txt, NoORFdataA1.txt,LMNmaxdataA1.txt, NoORFdataB1.txt, LMNmaxdataB1.txt, main.c, functions.c, functions.h, datain.c, datain.c, datain.h, print.c, print.h                                                                                                                    
                                                                                                                                                             
transfer_files = ALWAYS                                                                                                                                      
                                                                                                                                                             
Queue                                                                                                                                                        
" > something.classad

chmod 755 ${name}.sh
condor_submit something.classad
sleep 1

                            

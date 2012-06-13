#!/bin/bash           
                                                              
step1="`pwd`"
data="`basename $step1`"
step2="`dirname $step1`"
step3="`basename $step2`"
data=$step3"_"$data
model="SHM-t-df3_31MAY_THELASTdf5"

ORF="5000"

burn=$1
iter="10000"
thin="10"

name="CC_${model}_${data}_${ORF}_${burn}_${iter}_${thin}"
make
echo "#!/bin/sh                                                                                                                                              
hostname                                                                                                                                                                                                                                                                                          
                                                                              
date                                                                                                                                                         
chmod 755 main                                                                                                                                               
./main $burn $iter $thin $ORF >${name}.R                                                                                                                             
date                     " > ${name}.sh

echo "#####################                                                                                                                                  
                                                                                                                                                             
# My First Submission                                                                                                                                        
                                                                                                                                                             
#####################                                                                                                                                        
                                                                                                                                                             
Executable = ${name}.sh                                                                                                                              
                                                                                                                                                             
#Arguments = I bring you tidings of joy                                                                                                                      
                                                                                                                                                             
Universe = vanilla                                                                                                                                           
                                                                                                                                                             
#Requirements = OpSys == “LINUX” && Arch ==“X86_64”                                                                                                          
                                                                                                                                                             
Error = ${name}.err                                                                                                                                     
                                                                                                                                                             
Output = ${name}.out                                                                                                                                    
                                                                                                                                                             
Log = ${name}.log                                                                                                                                       
                                                                                                                                                             
transfer_input_files = main,headers.h, Makefile, priors.txt, ydata.txt, xdata.txt, NoTIMEdata.txt, NoORFdata.txt,LMNmaxdata.txt, main.c, functions.c, functions.h, datain.c, datain.h, print.c, print.h                                                                                                                               
                                                                                                                                                             
transfer_files = ALWAYS                                                                                                                                      
                                                                                                                                                             
Queue                                                                                                                                                        
" > something.classad

chmod 755 ${name}.sh
condor_submit something.classad
sleep 1

                            

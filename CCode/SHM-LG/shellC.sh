#!/bin/bash            
data="AdamCdc13-1" 
model="SHM"                                                                                                                                     
ORF="100"
name="CC_${model}_${data}_${ORF}"
burn="500000"
iter="100000"
thin="10"

cp ~/Model-in-C/CCode/Data/${model}/${data}/* ~/Model-in-C/CCode/${model}/

echo "#!/bin/sh                                                                                                                                              
hostname                                                                                                                                                     
make                                                                                                                                                         
sleep 10                                                                                                                                                     
date                                                                                                                                                         
chmod 755 main                                                                                                                                               
./main $burn $iter $thin $ORF >${name}.R                                                                                                                             
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
                                                                                                                                                             
transfer_input_files = main, Makefile, ydata.txt, xdata.txt, NoTIMEdata.txt, NoORFdata.txt,LMNmaxdata.txt, main.c, functions.c, functions.h, datain.c, datain.c datain.h, print.c, print.h                                                                                                                               
                                                                                                                                                             
transfer_files = ALWAYS                                                                                                                                      
                                                                                                                                                             
Queue                                                                                                                                                        
" > something.classad

chmod 755 ${name}.sh
condor_submit something.classad
sleep 1

                            

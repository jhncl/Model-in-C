#!/bin/bash  

step1="`pwd`"
data="`basename $step1`"
step2="`dirname $step1`"
step3="`basename $step2`"
data=$step3_$data
burn="100000"
iter="10000"
thin="10"

model="IHM" 
                                                                                             name="CC_${model}_${data}_${burn}_${iter}_${thin}"
make                                                                                                                                 
echo "#!/bin/sh                                                                                                                                              
hostname                                                                        

                                                                                                                                                                                                                                                                                                                                                                                     
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
                                                                                                                                                             
Error = ${name}.err                                                                                                                                     
                                                                                                                                                             
Output = ${name}.out                                               

Log = ${name}.log                                                                                                                                       
                                                                                                                                                             
transfer_input_files = priors.txt,main,headers.h, Makefile, dataA2.txt, dataB2.txt, NoORFdataA1.txt,LMNmaxdataA1.txt, NoORFdataB1.txt, LMNmaxdataB1.txt, main.c, functions.c, functions.h, datain.c, datain.c, datain.h, print.c, print.h                                                                                                                               
                                                                                                                                                             
transfer_files = ALWAYS                                                                                                                                      
                                                                                                                                                             
Queue                                                                                                                                                        
" > something.classad

chmod 755 ${name}.sh
condor_submit something.classad
sleep 1

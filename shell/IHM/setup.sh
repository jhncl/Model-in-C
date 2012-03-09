#!/bin/bash            
IFS=$'\n'
folder=`echo "$1_$3"`
mkdir $folder
folder_inner=`echo "$2_$4"`
mkdir $folder/$folder_inner
echo ".R file? $1 $2";
cd ~/QFADatasets/SHM/$1/$2
echo ls -lf CC_SHM_*.R
read inputline
whatA="$inputline"

echo ".R file? $3 $4";
cd ~/QFADatasets/SHM/$3/$4
echo ls -lf CC_SHM_*.R
read inputline
whatB="$inputline"
cd ~/QFADatasets/SHM/$1/$2/
echo "
aa<-read.table(\"$whatA\",header=T)                                     
dataA2=colMeans(aa)                                                        
write.table(dataA2,file=\"A2_$whatA\",row.names=F,col.names=F)                   " >colmeans.R
R CMD BATCH colmeans.R
mv ~/QFADatasets/SHM/$1/$2/A2_$whatA ~/QFADatasets/IHM/$folder/$folder_inner/dataA2.txt

echo "Half Way"

cd ~/QFADatasets/SHM/$3/$4/
echo "
aa<-read.table(\"$whatB\",header=T)                                     
dataB2=colMeans(aa)                                                        write.table(dataB2,file=\"B2_$whatB\",row.names=F,col.names=F)                   " >colmeans.R
R CMD BATCH colmeans.R
mv ~/QFADatasets/SHM/$3/$4/B2_$whatB ~/QFADatasets/IHM/$folder/$folder_inner/dataB2.txt
echo "Almost"

cp ~/QFADatasets/SHM/$1/$2/LMNmaxdata.txt ~/QFADatasets/IHM/$folder/$folder_inner/LMNmaxdataA1.txt

cp ~/QFADatasets/SHM/$1/$2/NoORFdata.txt ~/QFADatasets/IHM/$folder/$folder_inner/NoORFdataA\
1.txt
cp ~/QFADatasets/SHM/$3/$4/LMNmaxdata.txt ~/QFADatasets/IHM/$folder/$folder_inner/LMNmaxdataB1.txt

cp ~/QFADatasets/SHM/$3/$4/NoORFdata.txt ~/QFADatasets/IHM/$folder/$folder_inner/NoORFdataB1.txt
cp ~/QFADatasets/SHM/$1/M_SHM_FULL_$2.RData ~/QFADatasets/IHM/$folder/$folder_inner/M_SHM_FULL_$2_A1.RData
cp ~/QFADatasets/SHM/$3/M_SHM_FULL_$4.RData ~/QFADatasets/IHM/$folder/$folder_inner/M_SHM_FULL_$2_B1.RData

cp ~/QFADatasets/SHM/$3/$4/*.RData ~/QFADatasets/IHM/$folder/$folder_inner/.
cp ~/QFADatasets/SHM/$3/$4/*.RData ~/QFADatasets/IHM/$folder/$folder_inner/.

cp ~/QFADatasets/SHM/priors_IHM.txt ~/QFADatasets/IHM/$folder/$folder_inner/priors.txt
cp ~/2/Model-in-C/CCode/IHM/*.c ~/QFADatasets/IHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/IHM/*.h ~/QFADatasets/IHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/IHM/main ~/QFADatasets/IHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/IHM/Makefile ~/QFADatasets/IHM/$folder/$folder_inner/.
cp ~/2/Model-in-C/CCode/IHM/shellC.sh ~/QFADatasets/IHM/$folder/$folder_inner/.

cd ~/QFADatasets/IHM/$folder/$folder_inner/
make
echo "Done"
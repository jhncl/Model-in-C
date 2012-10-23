
library(qfaBayes,lib="~/R")
 library(qfa,lib="~/R")
library(rjags,lib="~/R")

a=read.delim("dataA1.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
b=read.delim("dataB1.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
TreatA=27
TreatB=27

filename=paste("M_JHM_FULL","_",TreatA,"_",TreatB,sep="")


qfa.variables(a)

Screen<-as.character(unique(a$Screen.Name))
MPlate<-as.character(unique(a$MasterPlate.Number))

a<-funcREMOVE(a,Screen,TreatA,MPlate)

a<-a[!a$Row==1,]
a<-a[!a$Row==16,]
a<-a[!a$Col==1,]
a<-a[!a$Col==24,]


Row<-a$Row
Col<-a$Col
for (i in 1:nrow(a)){
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

a$ID<-paste(a$Barcode,a$MasterPlate.Number,Row,Col,sep="")

ORFuni=unique(a$ORF)
########
#funcIDlist<-function(x){
#a$ID[a$ORF==x]
#}
#funcStrip<-function(x,i){x[1:i]}
#IDstrip=sapply(ORFuni,funcIDlist)
#IDstrip=sapply(IDstrip,unique)
#IDstrip=lapply(IDstrip,i=56,funcStrip)
#IDstrip=unlist(IDstrip)
#IDstrip=na.omit(IDstrip)
#a<-a[a$ID%in%IDstrip,]
#########
######################
if(file.exists("strip_list.txt")){
strip_list<-read.delim("strip_list.txt",header=T)
a<-a[!a$ORF%in%strip_list[,1],]
ORFuni<-unique(a$ORF)
}



a<-a[order(a$ORF,a$ID,a$Expt.Time), ]
ORFuni=unique(a$ORF)########

qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))
MPlate<-as.character(unique(b$MasterPlate.Number))
b<-funcREMOVE(b,Screen,TreatB,MPlate)

b<-b[!b$Row==1,]
b<-b[!b$Row==16,]
b<-b[!b$Col==1,]
b<-b[!b$Col==24,]


Row<-b$Row
Col<-b$Col
for (i in 1:nrow(b)){
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

b$ID<-paste(b$Barcode,b$MasterPlate.Number,Row,Col,sep="")
ORFuni_b<-unique(b$ORF)
####
sum(rep(1,length(ORFuni))[ORFuni==ORFuni_b])/length(ORFuni)
####

########
#funcIDlist<-function(x){
#b$ID[b$ORF==x]
#}
#funcStrip<-function(x,i){x[1:i]}
#IDstrip=sapply(ORFuni,funcIDlist)
#IDstrip=sapply(IDstrip,unique)
#IDstrip=lapply(IDstrip,i=56,funcStrip)
#IDstrip=unlist(IDstrip)
#IDstrip=na.omit(IDstrip)
#b<-b[b$ID%in%IDstrip,]
#########

if(file.exists("strip_list.txt")){
strip_list<-read.delim("strip_list.txt",header=T)
b<-b[!b$ORF%in%strip_list[,1],]
ORFuni_b<-unique(a$ORF)
}

#####################################
b<-b[order(b$ORF,b$ID,b$Expt.Time), ]
ORFuni_b<-unique(b$ORF)
sum(rep(1,length(ORFuni))[ORFuni==ORFuni_b])/length(ORFuni)


#a<-funcIDORDER(a)
IDuni<-unique(a$ID)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
if(sum(gene=="0")){
gene[gene=="0"]=ORFuni[gene=="0"]
}

N<-length(ORFuni);M=Ma=length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))

#b<-funcIDORDER(b)
IDuni<-unique(b$ID)
###



#gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);M=Mb=length(IDuni)#?
NoORF_b<-unlist(lapply(ORFuni,funcNoORF,data=b))#no of repeats each orf
NoTime_b<-c(0,unlist(lapply(IDuni,funcNoTime,data=b)))# 0+ no of time each repeat
NoSum_b<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_b)))


dimr<-max(NoORF_a,NoORF_b);dimc<-max(NoTime_a,NoTime_b)

funcXY_J<-function(data,data_b,Ma,Mb,N,NoTime_vec,NoSum_vec,NoTime_vec_b,NoSum_vec_b,dimr,dimc){
XY<-unlist(lapply(1:Ma,funcRowRep,NoTime_vec=NoTime_vec,data_vec=data,dimr,dimc))
XY<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec,data_vec=XY,dimr,dimc))
XY_b<-unlist(lapply(1:Mb,funcRowRep,NoTime_vec=NoTime_vec_b,data_vec=data_b,dimr,dimc))
XY_b<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec_b,data_vec=XY_b,dimr,dimc))
dim<-c(dimc,dimr,N,2)
XY<-funcARRAYTRANS_J(c(XY,XY_b),dim)
XY
}

y<-funcXY_J(a$Growth,b$Growth,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)
x<-funcXY_J(a$Expt.Time,b$Expt.Time,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)

QFA.I<-list("NoORF"=cbind(NoORF_a,NoORF_b),"NoTime"=cbind(NoTime_a,NoTime_b)[-1,],"NoSum"=cbind(NoSum_a,NoSum_b),"N"=N,"Ma"=Ma,"Mb"=Mb,"gene"=gene,SHIFT=c(0,max(NoSum_a,NoSum_b))
)
Scaling=FALSE######
if (Scaling==TRUE){y<-funcSCALING(rbind(a,b),y)}
QFA.D<-list(x=x,y=y)


x[is.na(x)]=-999
y[is.na(y)]=-999
xx<-aperm(x[,,,1],c(2,1,3))
yy<-aperm(y[,,,1],c(2,1,3))
write.table(file="xdataA1.txt",c(xx))
write.table(file="ydataA1.txt",c(yy))

write.table(file="NoORFdataA1.txt",c(NoORF_a))
write.table(file="NoTIMEdataA1.txt",c(NoTime_a)[-1])
write.table(file="LMNmaxdataA1.txt",c(N,max(NoORF_a),max(NoTime_a),length(y)/2,length(NoTime_a[-1])))

xx<-aperm(x[,,,2],c(2,1,3))
yy<-aperm(y[,,,2],c(2,1,3))
write.table(file="xdataB1.txt",c(xx))
write.table(file="ydataB1.txt",c(yy))

write.table(file="NoORFdataB1.txt",c(NoORF_b))
write.table(file="NoTIMEdataB1.txt",c(NoTime_b)[-1])
write.table(file="LMNmaxdataB1.txt",c(N,max(NoORF_b),max(NoTime_b),length(y)/2,length(NoTime_b[-1])))

save.image(paste(filename,".RData",sep=""))
stop()

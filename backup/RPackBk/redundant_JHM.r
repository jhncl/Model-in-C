data("Adam_cdc-1_SDLV2_REP1")
data("cdc13-1_rad9D_SDLv2_Rpt1")

TreatA=TreatB=27
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

ORFuni=unique(a$ORF)########
funcIDlist<-function(x){
a$ID[a$ORF==x]
}

funcStrip<-function(x,i){x[1:i]}
IDstrip=sapply(ORFuni,funcIDlist)
IDstrip=sapply(IDstrip,unique)
IDstrip=lapply(IDstrip,i=8,funcStrip)
IDstrip=unlist(IDstrip)
IDstrip=na.omit(IDstrip)
a<-a[a$ID%in%IDstrip,]
#########
######################


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

####

########
funcIDlist<-function(x){
b$ID[b$ORF==x]
}

funcStrip<-function(x,i){x[1:i]}
IDstrip=sapply(ORFuni_b,funcIDlist)
IDstrip=sapply(IDstrip,unique)
IDstrip=lapply(IDstrip,i=8,funcStrip)
IDstrip=unlist(IDstrip)
IDstrip=na.omit(IDstrip)
b<-b[b$ID%in%IDstrip,]
#########
#####################################
b<-b[order(b$ORF,b$ID,b$Expt.Time), ]
ORFuni_b<-unique(b$ORF)

sum(rep(1,length(ORFuni))[ORFuni==ORFuni_b])/length(ORFuni)

#a<-funcIDORDER(a)
IDuni<-unique(a$ID)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))

#b<-funcIDORDER(b)
IDuni<-unique(b$ID)
###



#gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);#?
NoORF_b<-unlist(lapply(ORFuni,funcNoORF,data=b))#no of repeats each orf
NoTime_b<-c(0,unlist(lapply(IDuni,funcNoTime,data=b)))# 0+ no of time each repeat
NoSum_b<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_b)))


dimr<-max(NoORF_a,NoORF_b);dimc<-max(NoTime_a,NoTime_b)



Ma<-length(unique(a$ID))
Mb<-length(unique(b$ID))
M<-max(Ma,Mb)

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

QFA.I<-list("NoORF"=cbind(NoORF_a,NoORF_b),"NoTime_a"=NoTime_a[-1],"NoTime_b"=NoTime_b[-1],"NoSum"=cbind(NoSum_a,NoSum_b),"N"=N,"M"=M,"gene"=gene,SHIFT=c(0,max(NoSum_a,NoSum_b))
)
Scaling=FALSE######
if (Scaling==TRUE){y<-funcSCALING(rbind(a,b),y)}
QFA.D<-list(x=x,y=y)

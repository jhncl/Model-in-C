data("Adam_cdc-1_SDLV2_REP1")
CAPN=6000
TREAT=unique(a$Treatment)
qfa.variables(a)
Treat=27
Screen<-unique(a$Screen.Name)
filename=paste("M_SHM_FULL","_",Treat,sep="")
MPlate<-(unique(a$MasterPlate.Number))
a<-funcREMOVE(a,Screen,Treat,MPlate)
a<-a[!a$Row==1,]
a<-a[!a$Row==16,]
a<-a[!a$Col==1,]
a<-a[!a$Col==24,]

Row<-paste(a$Row)
Col<-paste(a$Col)
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
IDstrip=lapply(IDstrip,unique)
IDstrip=lapply(IDstrip,i=8,funcStrip)
IDstrip=unlist(IDstrip)
IDstrip=na.omit(IDstrip)
a<-a[a$ID%in%IDstrip,]
#########

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]


Scaling=FALSE
IDuni<-unique(a$ID)
ORFuni=unique(a$ORF)########


gene<-unlist(lapply(ORFuni,funcGENE,data=a))
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))
dimr<-max(NoORF_a);dimc<-max(NoTime_a)
y<-funcXY(a$Growth,M,N,NoTime_a,NoSum_a,dimr,dimc)
x<-funcXY(a$Expt.Time,M,N,NoTime_a,NoSum_a,dimr,dimc)


QFA.I<-list("NoORF"=c(NoORF_a),"NoTime"=c(NoTime_a)[-1],"NoSum"=c(NoSum_a),"N"=N,"M"=M,"gene"=gene)
if (Scaling==TRUE){y<-funcSCALING(a,y)}
QFA.D<-list(y=y,x=x,ORFuni=ORFuni)

x[is.na(x)]=-999
y[is.na(y)]=-999
xx<-aperm(x,c(2,1,3))
yy<-aperm(y,c(2,1,3))

#################################################

 QFA.P<-list(

 sigma_K=7,               phi_K=0.1,
 eta_K_o=10,               psi_K_o=0.01,

 sigma_r=5,               phi_r=0.111,
 eta_r_o=6,               psi_r_o=0.111,

 eta_nu=-103.87,              psi_nu=0.1,

 K_mu=-17.95,     eta_K_p=0.25,
 r_mu=0.759,            eta_r_p=0.25,
 nu_mu=50,           eta_nu_p=0.0004,
 P_mu=-25.131,        eta_P=0.25
)


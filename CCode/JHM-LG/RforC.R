library(qfa,lib="~/R")

 library(qfaBayes,lib="~/R")



 library(qfaBayes)

Control<-c("Adam_cdc13-1_SDLV2_REP1.txt")
DescripControl<-"ExptDescriptionCDC13.txt"
a<-rod.read(files=Control[1],inoctimes=DescripControl)

Row<-paste(a$Row)
Col<-paste(a$Col)
for (i in 1:nrow(a)){
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}
a$ID<-paste(a$Barcode,a$MasterPlate.Number,Row,Col,sep="")
a<-a[order(a$ORF,a$ID,a$Expt.Time), ]


qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))
Treat<-27
MPlate<-15
a<-funcREMOVE(a,Screen,Treat,MPlate)


Scaling=FALSE
IDuni<-unique(a$ID)
ORFuni<-unique(a$ORF)


gene<-unlist(lapply(ORFuni,funcGENE,data=a))
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))
dimr<-max(NoORF_a);dimc<-max(NoTime_a)
y<-funcXY(a$Growth,M,N,NoTime_a,NoSum_a,dimr,dimc)
x<-funcXY(a$Expt.Time,M,N,NoTime_a,NoSum_a,dimr,dimc)
y<-y/max(y[!is.na(y)])###

QFA.I<-list("NoORF"=c(NoORF_a),"NoTime"=c(NoTime_a)[-1],"NoSum"=c(NoSum_a),"N"=
N,"M"=M,"gene"=gene)
if (Scaling==TRUE){y<-funcSCALING(a,y)}
QFA.D<-list(y=y,x=x,ORFuni=ORFuni)

x[is.na(x)]=-999
y[is.na(y)]=-999
xx<-aperm(x,c(2,1,3))
yy<-aperm(y,c(2,1,3))
write.table(file="xdata.txt",c(xx))
write.table(file="ydata.txt",c(yy))

write.table(file="NoORFdata.txt",c(NoORF_a))
write.table(file="NoTIMEdata.txt",c(NoTime_a)[-1])

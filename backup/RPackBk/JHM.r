#### Joint Hierachical Logistic Curve Model ####
qfa.Joint<-function(Control,Query,Scaling,iter,upd,thin,PlotOutput=TRUE,work,CustomModel=FALSE){
a<-Control
b<-Query
a<-funcIDORDER(a)
IDuni<-unique(a$ID)
ORFuni<-unique(a$ORF)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))

b<-funcIDORDER(b)
IDuni<-unique(b$ID)
ORFuni<-unique(b$ORF)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);M<-length(IDuni)#?
NoORF_b<-unlist(lapply(ORFuni,funcNoORF,data=b))#no of repeats each orf
NoTime_b<-c(0,unlist(lapply(IDuni,funcNoTime,data=b)))# 0+ no of time each repeat
NoSum_b<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_b)))


dimr<-max(NoORF_a,NoORF_b);dimc<-max(NoTime_a,NoTime_b)

Ma<-length(unique(a$ID))
Mb<-length(unique(b$ID))
M<-max(Ma,Mb)
y<-funcXY_J(a$Growth,b$Growth,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)
x<-funcXY_J(a$Expt.Time,b$Expt.Time,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)


QFA.I<-list("NoORF"=cbind(NoORF_a,NoORF_b),"NoTime"=cbind(NoTime_a,NoTime_b)[-1,],"NoSum"=cbind(NoSum_a,NoSum_b),"N"=N,"M"=M,"gene"=gene,SHIFT=c(0,max(NoSum_a,NoSum_b))
)

if (Scaling==TRUE){y<-funcSCALING(rbind(a,b),y)}
QFA.D<-list(x=x,y=y)
if (!(CustomModel==FALSE)){source(CustomModel)} else {funcMODELJoint()}
QFA.P<-funcPRIORS_J(CustomModel)

samp<-funcFITandUPDATE_J(QFA.I,QFA.D,QFA.P,iter,upd,thin)
QFA.O<-funcPosterior_J(samp,N,M,iter,thin,upd)

QFA<-c(QFA.O,QFA.I,QFA.D,QFA.P)
if(PlotOutput==TRUE){qfaplots.J(QFA,work)}
return(QFA)
}




### Joint Hierachical Logistic Curve Model Plots to Pdf###
qfaplots.J<-function(QFA,work,CustomInteractionDef=FALSE){
Treat="FIX"#######
samp<-QFA$samp
iter<-QFA$iter
thin<-QFA$thin

y<-QFA$y
x<-QFA$x

N<-QFA$N
M<-QFA$M
NoSum<-QFA$NoSum
NoORF<-QFA$NoORF
NoTime<-QFA$NoTime
gene<-QFA$gene
SHIFT<-QFA$SHIFT


K_s<-QFA$K_s
r_s<-QFA$r_s
PO_s<-QFA$PO_s
beta<-QFA$beta
tau_s<-QFA$tau_s
delta<-QFA$delta
alpha_ij_sd<-QFA$alpha_ij_sd
gamma_ij_sd<-QFA$gamma_ij_sd
alpha<-QFA$alpha
gamma<-QFA$gamma
alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij
p<-QFA$p
alpha_a<-QFA$alpha_a
alpha_b<-QFA$alpha_b
gam_b<-QFA$gam_b
omega_b<-QFA$omega_b



namesamp<-QFA$namesamp
K<-QFA$K
K_i<-QFA$K_i
K_ij<-QFA$K_ij
PO<-QFA$PO
k_tau<-QFA$k_tau
r<-QFA$r
r_i<-QFA$r_i
r_ij<-QFA$r_ij
r_tau<-QFA$r_tau
taui<-QFA$taui
tau<-QFA$tau
gam<-QFA$gam
omega<-QFA$omega
nu<-QFA$nu
nuc<-QFA$nuc
gamdelt<-QFA$gamdelt
omegadelt<-QFA$omegadelt
delta<-QFA$delta


A1<-QFA$alpha[1]
A2<-QFA$alpha[2]
B1<-QFA$bet[1]
B2<-QFA$bet[2]
sig<-sum(rep(1,N)[delta>0.5])
order<-order(1-delta)
vecorder<-order(1-delta)[1:sig]

################################################
print("Plots")
################################################

ylimmin<-0
ylimmax<-max(na.omit(as.numeric(y)))
xlimmin<-0
xlimmax<-max(na.omit(as.numeric(x)))

if(CustomInteractionDef==FALSE){
Mu<-funcInterDefMDRMDP(QFA)} else {
source(CustomInteractionDef)
Mu<-funcCustomInterDef(QFA)
}

Mu_a<-Mu[1:sum(NoORF[,1])]
veca<-vecb<-matrix(NA,N,max(NoORF))
Mu_b<-Mu[(1+sum(NoORF[,1])):sum(NoORF)]
mu_a=mu_b=0
for (i in 1:N){
mu_a[i]<-mean(Mu_a[(1+NoSum[i,1]):NoSum[i+1,1]])
mu_b[i]<-mean(Mu_b[(1+NoSum[i,2]):NoSum[i+1,2]])
}
limmin<-0
limmax<-max(na.omit(Mu))
###########################################
print("plot fitted with Conditioning on delta=1")
###########################################
pdf(paste("Plots_Inter",work,".pdf",sep=""))
plot(1,type="n",main=paste("Treatment",Treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control",ylab="Query",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=2)
abline(lm(c(colMeans(y[,,,2],na.rm=TRUE))~0+c(colMeans(y[,,,1],na.rm=TRUE))),col="grey",lty=3)
if (sum((1:N)[gene=="HIS3"])==1){
lines(c(mu_a[gene=="HIS3"],mu_b[gene=="HIS3"]),c(-1000,1000),lwd=2)
lines(c(-1000,1000),c(mean(mu_a[gene=="HIS3"]),mean(mu_b[gene=="HIS3"])),lwd=2)}
i=1:N
points(mu_a[i],mu_b[i],main=paste("Treatment",Treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single (=Alpha1*mu_i)",ylab="Double (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]>0]
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]<=0]  
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
legend(1,limmax, c("1-1","simple Lin Reg"), cex=0.5,col=c("grey","grey","black"), lty=c(1,2,3,1))
if (sum((1:N)[gene=="HIS3"])==1){
legend(1,limmax, c("1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("grey","grey","black"), lty=c(2,3,1))
}


plot(1,type="n",main=paste("Treatment",Treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control",ylab="Query",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=2)
abline(lm(c(colMeans(y[,,,1],na.rm=TRUE))~0+c(colMeans(y[,,,2],na.rm=TRUE))),col="grey",lty=3)
if (sum((1:N)[gene=="HIS3"])==1){
lines(c(mu_a[gene=="HIS3"],mu_b[gene=="HIS3"]),c(-1000,1000),lwd=2)
lines(c(-1000,1000),c(mean(mu_a[gene=="HIS3"]),mean(mu_b[gene=="HIS3"])),lwd=2)}
i=1:N
points(mu_a[i],mu_b[i],main=paste("Treatment",Treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single (=Alpha1*mu_i)",ylab="Double (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]>0]
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]<=0]  
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
legend(1,limmax, c("1-1","simple Lin Reg"), cex=0.5,col=c("grey","grey","black"), lty=c(2,3,1))
if (sum((1:N)[gene=="HIS3"])==1){
legend(1,limmax, c("1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("grey","grey","black"), lty=c(2,3,1))
}

i=1:N
plot(1,type="n",main=paste("Treatment",Treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control (=Alpha1*mu_i)",ylab="Query (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
lines(c(0,1000),c(0,1000),lwd=2)
points(mu_a[i],mu_b[i],main=paste("Treatment",Treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single (=Alpha1*mu_i)",ylab="Double (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]>0]
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]<=0]  
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=1:N
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
dev.off()

pdf(paste("Plots_M",work,".pdf",sep=""))
################################################
print("Master Curve")
################################################
plot(x,y,main="Master Curve",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(0,max(na.omit(c(x)))),ylim=c(ylimmin,ylimmax))
curve((K*PO*exp(r*x))/(K+PO*(exp(r*x)-1)), 0, max(na.omit(c(x))),add=TRUE,col=1)
################################################
print("ORF Curves")
################################################
plot(x[,,,1],y[,,,1],main="ORF Curves",xlab="Time (days)", ylab="Culture Domensity (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,,2],y[,,,2],main="ORF Curves",xlab="Time (days)", ylab="Culture Domensity (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
for (i in 1:N)
{
curve((K_i[i]*PO*exp(r_i[i]*x))/(K_i[i]+PO*(exp(r_i[i]*x)-1)), 0, 8,add=TRUE,col=2)
curve((A2*(K_i[i]+gamdelt[i])*PO*exp(B2*(r_i[i]+omegadelt[i])*x))/(A2*(K_i[i]+gamdelt[i])+PO*(exp(B2*(r_i[i]+omegadelt[i])*x)-1)), 0, 8,add=TRUE,col=3) 
}

################################################
print("Repeat Curves")
################################################
plot(x,y,main="Repeat Curves", xlab="Time (days)", ylab="Culture 
Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
for (i in 1:sum(NoSum[(N+1),]))
{
curve((K_ij[i]*PO*exp(r_ij[i]*x))/(K_ij[i]+PO*(exp(r_ij[i]*x)-1)), 0, 8,add=TRUE,col=i) 
}

dev.off()

###########################################
print("plots for individual Logistic curve fits")
###########################################
pdf(paste("Plots_M_indiv",work,".pdf",sep=""))
vecNoORF<-rbind(c(0,0),NoORF)
for (i in 1:N){
plot(x[,,i,1],y[,,i,1],col=2,main=paste(gene[i],"Repeat Curve"), xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i,2],y[,,i,2],col=3)
curve((K_i[i]*PO*exp(r_i[i]*x))/(K_i[i]+PO*(exp(r_i[i]*x)-1)), 0, 8,add=TRUE,col=2) 
curve((A2*(K_i[i]+gamdelt[i])*PO*exp(B2*(r_i[i]+omegadelt[i])*x))/(A2*(K_i[i]+gamdelt[i])+PO*(exp(B2*(r_i[i]+omegadelt[i])*x)-1)), 0, 8,add=TRUE,col=3) 

plot(x[,,i,1],y[,,i,1],col=2,main=paste(gene[i],"Repeat Curve"), xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i,2],y[,,i,2],col=3)
for (j in (1+sum(vecNoORF[1:i,1])):sum(vecNoORF[2:(i+1),1]))
{
curve((K_ij[j]*PO*exp(r_ij[j]*x))/(K_ij[j]+PO*(exp(r_ij[j]*x)-1)), 0, 8,add=TRUE,col=2) 
}
for (j in SHIFT[2]+((1+sum(vecNoORF[1:i,1])):sum(vecNoORF[2:(i+1),1])))
{
curve((K_ij[j]*PO*exp(r_ij[j]*x))/(K_ij[j]+PO*(exp(r_ij[j]*x)-1)), 0, 8,add=TRUE,col=3) 
}
}
dev.off()

pdf(paste("Plots_M_diag",work,".pdf",sep=""))
###########################################
print("Prior density")
###########################################
par(mfrow=c(4,2))
den<-matrix(0,2000,16)
den[,1]<-rgamma(2000,(K_s^2)/(alpha^2),K_s/(alpha^2))
den[,2]<-rgamma(2000,(K_s^2)/(alpha_i^2),K_s/(alpha_i^2))
den[,3]<-rgamma(2000,(K_s^2)/(alpha_ij^2),K_s/(alpha_ij^2))
den[,4]<-runif(2000,PO_s,beta)
den[,5]<-rgamma(2000,(r_s^2)/(gamma^2),r_s/(gamma^2))
den[,6]<-rgamma(2000,(r_s^2)/(gamma_i^2),r_s/(gamma_i^2))
den[,7]<-rgamma(2000,(r_s^2)/(gamma_ij^2),r_s/(gamma_ij^2))
den[,8]<-rgamma(2000,(tau_s^2)/(delta^2),tau_s/(delta^2))
den[,9]<-rgamma(2000,(alpha_ij^2)/(alpha_i^2),alpha_ij/(alpha_i^2))
for (i in 10:16){den[,i]<-rgamma(2000,(alpha_ij^2)/(alpha_i^2),alpha_ij/(alpha_i^2))
}
namesampden<-unique(substring(namesamp,1,4))
for (i in 1:ncol(den))
{
plot(density(den[,i]),paste(namesampden[i],"Prior Density"))
}
###########################################
print("Diagnostics trace acf density")
###########################################
par(mfrow=c(4,3))
for (i in c(1:(2*M+N+2),2*M+N+4,2*M+N+6,(2*M+4*N+7):(2*M+6*N+9),(2*M+7*N+18):(4*M+9*N+11)))
{
plot(as.numeric(samp[,i]),main=paste(namesamp[i],"Trace Top"),type="l")
plot(density(as.numeric(samp[,i])),main=paste(namesamp[i],"Density"))
t<-(1:ncol(den))[namesampden==substring(namesamp[i],1,4)]
lines(density(den[,t]),col=2)#####
acf(as.numeric(samp[,i]),main=paste(namesamp[i],"ACF"))
}
dev.off()
}

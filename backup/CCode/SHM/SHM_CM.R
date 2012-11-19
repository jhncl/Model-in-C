filename="M_SHM_A_FULL"

library(rjags)
library(qfa)
library(qfaBayes) 


Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")
DescripControl<-"ExptDescriptionCDC13.txt"

a<-rod.read(files=Control,inoctimes=DescripControl)


qfa.variables(a)
Screen<-unique(a$Screen.Name)
Treat<-27
MPlate<-(unique(a$MasterPlate.Number))
a<-funcREMOVE(a,Screen,Treat,MPlate)

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
IDstrip=lapply(IDstrip,i=6,funcStrip)
IDstrip=unlist(IDstrip)
IDstrip=na.omit(IDstrip)
a<-a[a$ID%in%IDstrip,]
#########

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]


Scaling=TRUE
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
write.table(file="LMNmaxdata.txt",c(N,max(NoORF_a),max(NoTime_a),length(y),length(NoTime_a[-1])))


#################################################
#################################################

 QFA.P<-list(

sigma_K=3,            phi_K=0.5, 
eta_K_o=0.5,          psi_K_o=1, 

 sigma_r=-1,               phi_r=0.1,
 eta_r_o=1,               psi_r_o=1,

 eta_nu=-1,              psi_nu=1,

 K_mu=log(0.2192928),     eta_K_p=1,
 r_mu=log(2.5),            eta_r_p=1,
 nu_mu=log(31),           eta_nu_p=1,
 P_mu=log(0.0002),        eta_P=1/0.01
)

write("
model {
      for (i in 1:N){
      	  for (j in 1:NoORF[i]){
	      	 for (l in 1:NoTime[(NoSum[i]+j)]){
			y[j,l,i] ~ dnorm(y.hat[j,l,i], exp(nu_l[i]))
			y.hat[j,l,i] <- (K_lm[(NoSum[i]+j)]*P*exp(r_lm[(NoSum[i]+j)]*x[j,l,i]))/(K_lm[(NoSum[i]+j)]+P*(exp(r_lm[(NoSum[i]+j)]*x[j,l,i])-1))
			}
		K_lm[(NoSum[i]+j)]<- exp(K_lm_L[(NoSum[i]+j)])
		K_lm_L[(NoSum[i]+j)] ~ dnorm(K_o_l[i],exp(tau_K_l[i]))
		r_lm[(NoSum[i]+j)]<- exp(min(3.5,r_lm_L[(NoSum[i]+j)]))
		r_lm_L[(NoSum[i]+j)] ~ dnorm(r_o_l[i],exp(tau_r_l[i]))
		}
	K_o_l[i] ~ dnorm( K_p, exp(sigma_K_o) )
	r_o_l[i] ~ dnorm( r_p, exp(sigma_r_o) )
	nu_l[i] ~ dnorm(nu_p,  exp(sigma_nu) )

	tau_K_l[i]~dnorm(sigma_K,phi_K)
	tau_r_l[i]~dnorm(sigma_r,phi_r)
	}

K_p ~ dnorm(K_mu,eta_K_p)
r_p ~ dnorm(r_mu,eta_r_p)
nu_p ~ dnorm(nu_mu,eta_nu_p)
P<-exp(P_L)
P_L ~ dnorm(P_mu,eta_P)

sigma_nu~dnorm(eta_nu,psi_nu)
sigma_K_o ~ dnorm(eta_K_o,psi_K_o)
sigma_r_o ~ dnorm(eta_r_o,psi_r_o)
}
","model1.bug")

 jags <- jags.model('model1.bug',
                    data = list('x' = QFA.D$x,
                                'y' = QFA.D$y,
                                'N' = QFA.I$N,
 'NoTime' = QFA.I$NoTime,
 'NoORF' = QFA.I$NoORF,
 'NoSum' = QFA.I$NoSum,
    'sigma_K'=QFA.P$sigma_K,                'phi_K'=QFA.P$phi_K,
    'sigma_r'=QFA.P$sigma_r,                'phi_r'=QFA.P$phi_r,
    'eta_K_o'=QFA.P$eta_K_o,                'psi_K_o'=QFA.P$psi_K_o,
    'eta_r_o'=QFA.P$eta_r_o,                'psi_r_o'=QFA.P$psi_r_o,
    'eta_nu'=QFA.P$eta_nu,               'psi_nu'=QFA.P$psi_nu,
    'K_mu'=QFA.P$K_mu,                   'eta_K_p'=QFA.P$eta_K_p,
    'r_mu'=QFA.P$r_mu,                   'eta_r_p'=QFA.P$eta_r_p,
    'nu_mu'=QFA.P$nu_mu,                  'eta_nu_p'=QFA.P$eta_nu_p,
    'P_mu'=QFA.P$P_mu,                   'eta_P'=QFA.P$eta_P
 ),
                    n.chains = 1,
                    n.adapt = 100)

samp<-coda.samples(jags,
 c('K_lm_L',            'tau_K_l',
    'K_o_l',            'sigma_K_o',  
    'K_p',
    'P_L',
    'r_lm_L',            'tau_r_l',
    'r_o_l',            'sigma_r_o',
    'r_p',
    'nu_l',             'sigma_nu',
    'nu_p'),
              1000,thin=1)
save(samp,file=paste(filename,"_F0.R",sep=""))


 update(jags,10000)
date()
samp<-coda.samples(jags,
 c('K_lm_L',            'tau_K_l',
    'K_o_l',            'sigma_K_o',  
    'K_p',
    'P_L',
    'r_lm_L',            'tau_r_l',
    'r_o_l',            'sigma_r_o',
    'r_p',
    'nu_l',             'sigma_nu',
    'nu_p'),
              20000,thin=20)

save(samp,file=paste(filename,"_F1.R",sep=""))


samp<-coda.samples(jags,
 c('K_lm_L',            'tau_K_l',
    'K_o_l',            'sigma_K_o',  
    'K_p',
    'P_L',
    'r_lm_L',            'tau_r_l',
    'r_o_l',            'sigma_r_o',
    'r_p',
    'nu_l',             'sigma_nu',
    'nu_p'),
              100000,thin=10)

save(samp,file=paste(filename,"_F2.R",sep=""))

stop()

L=N=4294
M=17308
pdf(file="SHMtestplot_C_4294.pdf")
#K_lm[%i]
for (i in 1:M){
j=i
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
}
dev.off()
#tau_K_l[%i]
j=M+1
for (i in (2*M+3*N+8):(2*M+4*N+7)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
}

#"K_o_l[%i] 
j=M+N+1
for (i in (M+1):(M+N)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
}

#sigma_K_o ");
i=2*M+3*N+5
j=M+2*N+1
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#K_p ");
i=M+1+N
j=M+2*N+2
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#"P_l ");
i=(M+N+2)
j=M+2*N+3
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#r_lm[%i] 
j=M+2*N+4
for (i in (M+2*N+4):(2*M+2*N+3)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#tau_r_l[%i] ",l);
j=2*M+2*N+4
for (i in (2*M+4*N+8):(2*M+5*N+7)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#r_o_l[%i] ",l);
j=2*M+3*N+4
for (i in (2*M+2*N+4):(2*M+3*N+3)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))j=j+1
}

#sigma_r_o ");
i=2*M+3*N+7
j=2*M+4*N+4
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#r_p ");
i=2*M+3*N+4
j=2*M+4*N+5
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#"nu_l[%i] ",l);
j=2*M+4*N+6
for (i in (M+N+3):(M+2*N+2)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_nu ");
i=2*M+3*N+6
j=2*M+5*N+6
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#nu_p ");
i=M+2*N+3
j=2*M+5*N+7
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
dev.off()

i=M+N+3;j=i+N+1;
colnames(samp)[i];names(aa)[j]




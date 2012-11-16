funcMCurveVar_LG<-function(QFA){
x<-QFA$x
samp<-QFA$samp
M<-QFA$M
N<-QFA$N
MQQU=MQQD=NA
vec<-seq(0,ceiling(max(x,na.rm=TRUE)),0.1)
QQ<-matrix(numeric(0),ncol=length(vec),nrow=nrow(samp))

for (l in 1:nrow(samp)){
PO=samp[l,(M+2*N+3)]
KK=exp(rnorm(1,(samp[l,1]),1/samp[l,N+2]^0.5))
rr=exp(rnorm(1,(samp[l,(M+2*N+5)]),1/samp[l,(M+3*N+6)]^0.5))
QQ[l,]<-(KK*PO*exp(rr*vec))/(KK+PO*(exp(rr*vec)-1))
}
for (j in 1:length(vec) ){
MQQU[j]<-quantile(QQ[,j],na.rm=TRUE)[4]
MQQD[j]<-quantile(QQ[,j],na.rm=TRUE)[2]
}
list(MQQU=MQQU,MQQD=MQQD,Time=vec)
}

funcICurveVar_LG<-function(QFA){
x<-QFA$x
y<-QFA$y

gene<-QFA$gene
ylimmin<-0
ylimmax<-max(na.omit(as.numeric(y)))
xlimmin<-0
xlimmax<-max(na.omit(as.numeric(x)))

samp<-QFA$samp
PO<-QFA$PO
M<-QFA$M
N<-QFA$N
MQQU=MQQD=NA
NoSum<-QFA$NoSum
NoORF<-QFA$NoORF
NoTime<-QFA$NoTime
K_i<-exp(QFA$K_i)
K_ij<-QFA$K_ij
r_i<-exp(QFA$r_i)
r_ij<-QFA$r_ij

for (i in 1:N){

 plot(-1,-1,main=paste(gene[i],"Curve",i),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i],y[,,i])
vec<-seq(0,ceiling(max(x,na.rm=TRUE)),0.1)
QQ<-matrix(numeric(0),ncol=length(vec),nrow=nrow(samp))
 for (l in 1:nrow(samp)){
  PO=samp[l,(M+2*N+3)]
  A=1/samp[l,M+N+3+i-1]^0.5
  B=samp[l,2*M+3*N+7+i-1]
  KK=exp(rnorm(1,samp[l,(2+i-1)],A))
  rr=exp(rnorm(1,samp[l,(M+2*N+6+i-1)],B))
QQ[l,]<-(KK*PO*exp(rr*vec))/(KK+PO*(exp(rr*vec)-1))
 }

for (j in 1:length(vec)){
MQQU[j]<-quantile(QQ[,j],na.rm=TRUE)[4]
MQQD[j]<-quantile(QQ[,j],na.rm=TRUE)[2]
}
MCurveVar<-list(MQQU=MQQU,MQQD=MQQD)
KK=K_i[i]
rr=r_i[i]
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)),xlimmin, xlimmax,add=TRUE) 
lines(vec,MCurveVar$MQQU,col=2)
lines(vec,MCurveVar$MQQD,col=2)

plot(-1,-1,main=paste(gene[i],"Repeat Curves"),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i],y[,,i])
for (j in 1:NoORF[i])
{
KK=K_ij[(j+NoSum[i])]
rr=r_ij[(j+NoSum[i])]
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=1+j) 
}
}
}





















funcModelVarPost_LG<-function(QFA,i){
K<-QFA$K
alpha<-QFA$alpha
K_i<-QFA$K_i
k_tau<-QFA$k_tau
K_ij<-QFA$K_ij
r<-QFA$r
gamma<-QFA$gamma
r_i<-QFA$r_i
gamma_i<-QFA$gamma_i
r_ij<-QFA$r_ij
NoSum<-QFA$NoSum

len<-length(K_ij[((1+NoSum[i]):NoSum[i+1])])
list(
A=rnorm(3*len,K,alpha),
B=rnorm(2*len,K_i[i],1/k_tau[i]^0.5),
C=log(K_ij[((1+NoSum[i]):NoSum[i+1])]),
D=rnorm(3*len,r,gamma),
E=rnorm(2*len,r_i[i],gamma_i),
F=log(r_ij[((1+NoSum[i]):NoSum[i+1])])
)
}

funcDen_LG<-function(sampsize,QFA){
K_s<-QFA$K_s
r_s<-QFA$r_s
PO_s<-QFA$PO_s
beta<-QFA$beta
tau_s<-QFA$tau_s
delta<-QFA$delta
alpha_i_tau<-QFA$alpha_i_tau
gamma_i_tau<-QFA$gamma_i_tau
alpha<-QFA$alpha
gamma<-QFA$gamma
alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij
alpha_ij_tau<-QFA$alpha_ij_tau
gamma_ij_tau<-QFA$gamma_ij_tau
delta_mu<-QFA$delta_mu

delta_sd<-QFA$delta_sd

den<-matrix(0,sampsize,14)
den[,1]<-(rnorm(sampsize,K_s,alpha))
den[,2]<-(rnorm(sampsize,K_s,alpha))
den[,3]<-(rnorm(sampsize,alpha_i,1/alpha_i_tau^0.5))
den[,4]<-exp(rnorm(sampsize,K_s,alpha))
den[,5]<-exp(rnorm(sampsize,alpha_ij,1/alpha_ij_tau^0.5))
den[,6]<-exp(rnorm(sampsize,PO_s,beta))
den[,7]<-exp(rnorm(sampsize,delta_mu,1/delta_sd^0.5))
den[,8]<-(rnorm(sampsize,r_s,gamma))
den[,9]<-(rnorm(sampsize,r_s,gamma))
den[,10]<-exp(rnorm(sampsize,gamma_i,1/gamma_i_tau^0.5))
den[,11]<-exp(rnorm(sampsize,r_s,gamma))
den[,12]<-exp(rnorm(sampsize,gamma_ij,1/gamma_ij_tau^0.5))
den[,13]<-rnorm(sampsize,tau_s,1/delta^0.5)
den[,14]<-exp(rnorm(sampsize,tau_s,1/delta^0.5))
den
}

funcPostPred_LG<-function(iter,QFA){

N<-QFA$N
M<-QFA$M
NoSum<-QFA$NoSum
NoORF<-QFA$NoORF
NoTime<-QFA$NoTime

K_s<-QFA$K_s
r_s<-QFA$r_s
PO_s<-QFA$PO_s
beta<-QFA$beta
tau_s<-QFA$tau_s
delta<-QFA$delta
alpha_i_tau<-QFA$alpha_i_tau
gamma_i_tau<-QFA$gamma_i_tau
alpha<-QFA$alpha
gamma<-QFA$gamma
alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij
alpha_ij_tau<-QFA$alpha_ij_tau
gamma_ij_tau<-QFA$gamma_ij_tau

namesamp<-QFA$namesamp
K<-QFA$K
K_i<-QFA$K_i
K_ij<-QFA$K_ij
PO<-QFA$PO
K_ij_tau<-QFA$K_ij_tau
r<-QFA$r
r_i<-QFA$r_i
r_ij<-QFA$r_ij
r_ij_tau<-QFA$r_ij_tau
taui<-QFA$taui
tau<-QFA$tau
delta_mu<-QFA$delta_mu

delta_sd<-QFA$delta_sd
#exp Kij tauij rij kvar i ij rvar i ij  po  

samp<-QFA$samp



postpred<-matrix(0,iter,length(namesamp))

postpred[,1]<-rnorm(iter,K_s,1/alpha^0.5)

for (i in 1:length(K_i)){
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c( Q,rnorm(1,(samp[l,2+i-1]),1/samp[l,N+2]^0.5))
}
A=mean(Q);B=sd(Q)
postpred[,2+i-1]<-(rnorm(iter,(A),B))
}

postpred[,(N+2)]<-exp(rnorm(iter,alpha_i,1/alpha_i_tau^0.5))

t=1
for (i in (N+3):(M+N+2)){
if(((i-(N+3))==c(NoSum[-1])[t])) t=t+1
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rnorm(1,(samp[l,2+t-1]),1/samp[l,M+N+3+t-1]^0.5))
}
A=mean(Q);B=sd(Q)
postpred[,i]<-exp(rnorm(iter,(A),B))
}

for (i in 1:length(K_ij_tau)){
j=(M+N+3)+i-1
postpred[,j]<-exp(rnorm(iter,(alpha_ij),1/alpha_ij_tau^0.5))
}

postpred[,(M+2*N+3)]<-exp(rnorm(iter,(PO_s),1/beta^0.5))

postpred[,(M+2*N+4)]<-exp(rnorm(iter,(delta_mu),1/delta_sd^0.5))

##
postpred[,(M+2*N+5)]<-rnorm(iter,(r_s),1/gamma^0.5)

for (i in 1:length(r_i)){
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rnorm(1,(samp[l,M+2*N+6+i-1]),1/samp[l,M+3*N+6]^0.5))
}
A=mean(Q);B=sd(Q)
postpred[,M+2*N+6+i-1]<-rnorm(iter,(A),B)
}

postpred[,M+3*N+6]<-exp(rnorm(iter,(gamma_i),1/gamma_i_tau^0.5))

t=1
for (i in (M+3*N+7):(2*M+3*N+6)){
if(((i-(M+3*N+7))==c(NoSum[-1])[t])) t=t+1
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rnorm(1,(samp[l,M+2*N+6+t-1]^2),samp[l,2*M+3*N+7+t-1]))
}
A=mean(Q);B=sd(Q)
postpred[,i]<-exp(rnorm(iter,(A),B))
}

for (i in 1:length(r_ij_tau)){
j=(2*M+3*N+7)+i-1
postpred[,j]<-exp(rnorm(iter,(gamma_ij),1/gamma_ij_tau^0.5))
}
##

postpred[,(2*M+4*N+7)]<-exp(rnorm(iter,(tau_s),1/delta^0.5))

Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rnorm(iter,(samp[l,(2*M+4*N+7)]),1/samp[l,(M+2*N+4)]^0.5))
}
A=mean(Q);B=sd(Q)
postpred[,(2*M+4*N+8):(2*M+5*N+7)]<-exp(rnorm(iter,(A),B))
postpred
}





funcPriorPlot_LG<-function(QFA.P){
K_s<-QFA.P$K_s
r_s<-QFA.P$r_s
PO_s<-QFA.P$PO_s
beta<-QFA.P$beta
tau_s<-QFA.P$tau_s
delta<-QFA.P$delta

alpha<-QFA.P$alpha
gamma<-QFA.P$gamma

alpha_i<-QFA.P$alpha_i
gamma_i<-QFA.P$gamma_i
alpha_i_tau<-QFA.P$alpha_i_tau
gamma_i_tau<-QFA.P$gamma_i_tau

alpha_ij<-QFA.P$alpha_ij
gamma_ij<-QFA.P$gamma_ij
alpha_ij_tau<-QFA.P$alpha_ij_tau
gamma_ij_tau<-QFA.P$gamma_ij_tau

delta_mu<-QFA.P$delta_mu
delta_sd<-QFA.P$delta_sd

sampsize=2000
pdf(paste("Plots_M_Priors_",".pdf",sep=""))

name<-c("K Prior","K_i_tau Prior","K_ij_tau Prior","r Prior","r_i_tau Prior","r_ij_tau Prior","PO Prior","tau Prior","delta_tau")
vec<-matrix(numeric(0),sampsize,9)
vec[,1]<-exp(rnorm(sampsize,(K_s),alpha))
vec[,2]<-exp(rnorm(sampsize,(alpha_i),alpha_i_tau))
vec[,3]<-exp(rnorm(sampsize,(alpha_ij),alpha_ij_tau))
vec[,4]<-exp(rnorm(sampsize,(r_s),gamma))
vec[,5]<-exp(rnorm(sampsize,(gamma_i),gamma_i_tau))
vec[,6]<-exp(rnorm(sampsize,(gamma_ij),gamma_ij_tau))
vec[,7]<-exp(rnorm(sampsize,(PO_s),beta))
vec[,8]<-exp(rnorm(sampsize,(tau_s),1/delta^0.5))
vec[,9]<-exp(rnorm(sampsize,(delta_mu),delta_sd))

par(mfrow=c(3,3))
for (i in 1:length(name)){
plot(density(vec[,i]),main=name[i])
}
par(mfrow=c(1,1))
for (i in 1:length(name)){
plot(density(vec[,i]),main=name[i])
mtext( paste(round(quantile(vec[,i])[1],digit=3),
round(quantile(vec[,i])[2],digits=3),
round(quantile(vec[,i])[3],digits=3),
round(quantile(vec[,i])[4],digits=3),
round(quantile(vec[,i])[5],digits=3)
))
}

#sd
name<-c("K Prior","K_i_tau Prior","K_ij_tau Prior","r Prior","r_i_tau Prior","r_ij_tau Prior","PO Prior","tau Prior","delta_tau")
vec<-matrix(numeric(0),sampsize,9)
vec[,1]<-exp(rnorm(sampsize,(K_s),alpha))
vec[,2]<-1/exp(rnorm(sampsize,(alpha_i),alpha_i_tau))^0.5
vec[,3]<-1/exp(rnorm(sampsize,(alpha_ij),alpha_ij_tau))^0.5
vec[,4]<-exp(rnorm(sampsize,(r_s),gamma))
vec[,5]<-1/exp(rnorm(sampsize,(gamma_i),gamma_i_tau))^0.5
vec[,6]<-exp(rnorm(sampsize,(gamma_ij),gamma_ij_tau))
vec[,7]<-exp(rnorm(sampsize,(PO_s),beta))
vec[,8]<-1/exp(rnorm(sampsize,(tau_s),1/delta^0.5))^0.5
vec[,9]<-1/exp(rnorm(sampsize,(delta_mu),delta_sd))^0.5

par(mfrow=c(3,3))
for (i in 1:length(name)){
plot(density(vec[,i]),main=name[i])
}
par(mfrow=c(1,1))
for (i in 1:length(name)){
plot(density(vec[,i]),main=name[i])
mtext( paste(round(quantile(vec[,i])[1],digits=3),
round(quantile(vec[,i])[2],digits=3),
round(quantile(vec[,i])[3],digits=3),
round(quantile(vec[,i])[4],digits=3),
round(quantile(vec[,i])[5],digits=3)
))
}

dev.off()

}



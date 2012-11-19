funcMCurveVar<-function(QFA){
x<-QFA$x
samp<-QFA$samp
M<-QFA$M
N<-QFA$N
MQQU=MQQD=NA
vec<-seq(0,ceiling(max(x,na.rm=TRUE)),0.1)
QQ<-matrix(numeric(0),ncol=length(vec),nrow=nrow(samp))
for (l in 1:nrow(samp)){
PO=samp[l,(M+2*N+3)]
KK=rgamma(1,(samp[l,1]^2)*samp[l,N+2],samp[l,1]*samp[l,N+2])
rr=rgamma(1,(samp[l,(M+2*N+5)]^2)*(samp[l,(M+3*N+6)]),samp[l,(M+2*N+5)]*(samp[l,(M+3*N+6)]))
QQ[l,]<-(KK*PO*exp(rr*vec))/(KK+PO*(exp(rr*vec)-1))
}

for (j in 1:length(vec) ){
MQQU[j]<-quantile(QQ[,j],na.rm=TRUE)[4]
MQQD[j]<-quantile(QQ[,j],na.rm=TRUE)[2]
}
list(MQQU=MQQU,MQQD=MQQD,Time=vec)
}

funcICurveVar<-function(QFA){
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
NoSum<-QFA$NoSum
NoORF<-QFA$NoORF
NoTime<-QFA$NoTime
MQQU=MQQD=NA
K_i<-QFA$K_i
K_ij<-QFA$K_ij
r_i<-QFA$r_i
r_ij<-QFA$r_ij

for (i in 1:N){

 plot(-1,-1,main=paste(gene[i],"Curve",i),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i],y[,,i])
vec<-seq(0,ceiling(max(x,na.rm=TRUE)),0.1)
QQ<-matrix(numeric(0),ncol=length(vec),nrow=nrow(samp))
 for (l in 1:nrow(samp)){
  PO=samp[l,(M+2*N+3)]
  A=samp[l,M+N+3+i-1]
  B=1/samp[l,2*M+3*N+7+i-1]^2
  KK=rgamma(1,(samp[l,(2+i-1)]^2)*A,samp[l,(2+i-1)]*A)
  rr=rgamma(1,(samp[l,(M+2*N+6+i-1)]^2)*B,samp[l,(M+2*N+6+i-1)]*B)
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


funcDen<-function(sampsize,QFA){
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
alpha_ij_tau<-QFA$alpha_ij
gamma_ij_tau<-QFA$gamma_ij_tau
delta_mu<-QFA$delta_mu

delta_sd<-QFA$delta_sd

den<-matrix(0,sampsize,14)
den[,1]<-rgamma(sampsize,(K_s^2)*alpha,K_s*alpha)
den[,2]<-rgamma(sampsize,(K_s^2)*alpha,K_s*alpha)
den[,3]<-rgamma(sampsize,(alpha_i^2)*alpha_i_tau,alpha_i*alpha_i_tau)
den[,4]<-rgamma(sampsize,(K_s^2)*alpha,K_s*alpha)
den[,5]<-rgamma(sampsize,(alpha_ij^2)*alpha_ij_tau,alpha_ij*alpha_ij_tau)
den[,6]<-rgamma(sampsize,(PO_s^2)*beta,PO_s*beta)
den[,7]<-rgamma(sampsize,(delta_mu^2)*delta_sd,delta_mu*delta_sd)
den[,8]<-rgamma(sampsize,(r_s^2)*gamma,r_s*gamma)
den[,9]<-rgamma(sampsize,(r_s^2)*gamma,r_s*gamma)
den[,10]<-rgamma(sampsize,(gamma_i^2)*gamma_i_tau,gamma_i*gamma_i_tau)
den[,11]<-rgamma(sampsize,(r_s^2)*gamma,r_s*gamma)
den[,12]<-rgamma(sampsize,(gamma_ij^2)*gamma_ij_tau,gamma_ij*gamma_ij_tau)
den[,13]<-rgamma(sampsize,(tau_s^2)*delta,tau_s*delta)
den[,14]<-rgamma(sampsize,(tau_s^2)*delta,tau_s*delta)
den
}

funcPostPred<-function(iter,QFA){

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
alpha_ij_tau<-QFA$alpha_ij_tau
gamma_ij_tau<-QFA$gamma_ij_tau
alpha_i_tau<-QFA$alpha_i_tau
gamma_i_tau<-QFA$gamma_i_tau
alpha<-QFA$alpha
gamma<-QFA$gamma
alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij

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

samp<-QFA$samp


postpred<-matrix(0,iter,length(namesamp))

postpred[,1]<-rgamma(iter,(K_s^2)*alpha,K_s*alpha)


for (i in 1:length(K_i)){
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rgamma(1,(samp[l,2+i-1]^2)*samp[l,N+2],samp[l,2+i-1]*samp[l,N+2]))
}
A=mean(Q);B=1/var(Q)
postpred[,2+i-1]<-rgamma(iter,(A^2)*B,A*B)
}

postpred[,(N+2)]<-rgamma(iter,(alpha_i^2)*alpha_i_tau,alpha_i*alpha_i_tau)

t=1
for (i in (N+3):(M+N+2)){
if(((i-(N+3))==c(NoSum[-1])[t])) t=t+1
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rgamma(1,(samp[l,2+t-1]^2)*samp[l,M+N+3+t-1],samp[l,2+t-1]*samp[l,M+N+3+t-1]))
}
A=mean(Q);B=1/var(Q)
postpred[,i]<-rgamma(iter,(A^2)*B,A*B)
}

for (i in 1:length(K_ij_tau)){
j=(M+N+3)+i-1
postpred[,j]<-rgamma(iter,(alpha_ij^2)*alpha_ij_tau,alpha_ij*alpha_ij_tau)
}

postpred[,(M+2*N+3)]<-rgamma(iter,(PO_s^2)*beta,PO_s*beta)

postpred[,(M+2*N+4)]<-rgamma(iter,(delta_mu^2)*delta_sd,delta_mu*delta_sd)

##
postpred[,(M+2*N+5)]<-rgamma(iter,(r_s^2)*gamma,r_s*gamma)

for (i in 1:length(r_i)){
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rgamma(1,(samp[l,M+2*N+6+i-1]^2)*samp[l,M+3*N+6],samp[l,M+2*N+6+i-1]*samp[l,M+3*N+6]))
}
A=mean(Q);B=1/var(Q)
postpred[,M+2*N+6+i-1]<-rgamma(iter,(A^2)*B,A*B)
}

postpred[,M+3*N+6]<-rgamma(iter,(gamma_i^2)*gamma_i_tau,gamma_i*gamma_i_tau)

t=1
for (i in (M+3*N+7):(2*M+3*N+6)){
if(((i-(M+3*N+7))==c(NoSum[-1])[t])) t=t+1
Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rgamma(1,(samp[l,M+2*N+6+t-1]^2)/samp[l,2*M+3*N+7+t-1]^2,samp[l,M+2*N+6+t-1]/samp[l,2*M+3*N+7+t-1]^2))
}
A=mean(Q);B=1/var(Q)
postpred[,i]<-rgamma(iter,(A^2)*B,A*B)
}

for (i in 1:length(r_ij_tau)){
j=(2*M+3*N+7)+i-1
postpred[,j]<-rgamma(iter,(gamma_ij^2)*gamma_ij_tau,gamma_ij*gamma_ij_tau)
}
##

postpred[,(2*M+4*N+7)]<-rgamma(iter,(tau_s^2)*delta,tau_s*delta)

Q<-numeric(0)
for (l in 1:nrow(samp)){
Q<-c(Q,rgamma(iter,(samp[l,(2*M+4*N+7)]^2)*samp[l,(M+2*N+4)],samp[l,(2*M+4*N+7)]*samp[l,(M+2*N+4)]))
}
A=mean(Q);B=1/var(Q)
postpred[,(2*M+4*N+8):(2*M+5*N+7)]<-rgamma(iter,(A^2)*B,A*B)
postpred
}

funcPriorPlot<-function(QFA.P,work){

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

#para
sampsize=2000
pdf(paste("Plots_M_Priors_",".pdf",sep=""))

name<-c("K Prior","K_i_tau Prior","K_ij_tau Prior","r Prior","r_i_tau Prior","r_ij_tau Prior","PO Prior","tau Prior","delta_tau")
vec<-matrix(numeric(0),sampsize,9)
vec[,1]<-rgamma(sampsize,(K_s^2)*alpha,K_s*alpha)
vec[,2]<-rgamma(sampsize,(alpha_i^2)*alpha_i_tau,alpha_i*alpha_i_tau)
vec[,3]<-rgamma(sampsize,(alpha_ij^2)*alpha_ij_tau,alpha_ij*alpha_ij_tau)
vec[,4]<-rgamma(sampsize,(r_s^2)*gamma,r_s*gamma)
vec[,5]<-rgamma(sampsize,(gamma_i^2)*gamma_i_tau,gamma_i*gamma_i_tau)
vec[,6]<-rgamma(sampsize,(gamma_ij^2)/gamma_ij_tau,gamma_ij/gamma_ij_tau)
vec[,7]<-rgamma(sampsize,(PO_s^2)*beta,PO_s*beta)
vec[,8]<-rgamma(sampsize,(tau_s^2)*delta,tau_s*delta)
vec[,9]<-rgamma(sampsize,(delta_mu^2)*delta_sd,delta_mu*delta_sd)

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
vec[,1]<-rgamma(sampsize,(K_s^2)*alpha,K_s*alpha)
vec[,2]<-1/rgamma(sampsize,(alpha_i^2)*alpha_i_tau,alpha_i*alpha_i_tau)^0.5
vec[,3]<-1/rgamma(sampsize,(alpha_ij^2)*alpha_ij_tau,alpha_ij*alpha_ij_tau)^0.5
vec[,4]<-rgamma(sampsize,(r_s^2)*gamma,r_s*gamma)
vec[,5]<-1/rgamma(sampsize,(gamma_i^2)*gamma_i_tau,gamma_i*gamma_i_tau)^0.5

vec[,6]<-rgamma(sampsize,(gamma_ij^2)/gamma_ij_tau,gamma_ij/gamma_ij_tau)

vec[,7]<-rgamma(sampsize,(PO_s^2)*beta,PO_s*beta)
vec[,8]<-1/rgamma(sampsize,(tau_s^2)*delta,tau_s*delta)^0.5
vec[,9]<-1/rgamma(sampsize,(delta_mu^2)*delta_sd,delta_mu*delta_sd)^0.5

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

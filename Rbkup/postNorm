funcDen_LG<-function(sampsize,QFA){
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

den<-matrix(0,sampsize,11)
den[,1]<-rnorm(sampsize,K_s,alpha)
den[,2]<-rnorm(sampsize,K_s,alpha)
den[,3]<-exp(rnorm(sampsize,K_s,alpha))
den[,4]<-runif(sampsize,PO_s,beta)
den[,5]<-rgamma(sampsize,(alpha_ij^2)/(alpha_ij_sd^2),alpha_ij/(alpha_ij_sd^2))
den[,6]<-rnorm(sampsize,r_s,gamma)
den[,7]<-rnorm(sampsize,r_s,gamma)
den[,8]<-exp(rnorm(sampsize,r_s,gamma))
den[,9]<-rgamma(sampsize,(gamma_ij^2)/(gamma_ij_sd^2),gamma_ij/(gamma_ij_sd^2))
den[,10]<-rgamma(sampsize,(tau_s^2)/(delta^2),tau_s/(delta^2))
den[,11]<-rgamma(sampsize,(tau_s^2)/(delta^2),tau_s/(delta^2))
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
alpha_ij_sd<-QFA$alpha_ij_sd
gamma_ij_sd<-QFA$gamma_ij_sd
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
k_tau<-QFA$k_tau
r<-QFA$r
r_i<-QFA$r_i
r_ij<-QFA$r_ij
r_tau<-QFA$r_tau
taui<-QFA$taui
tau<-QFA$tau





postpred<-matrix(0,iter,length(namesamp))

postpred[,1]<-rnorm(iter,K_s,alpha)

for (i in 1:length(K_i)){
j=(2)+i-1
postpred[,j]<-rnorm(iter,K,alpha_i)
}

t=1
for (i in (N+2):(M+N+1) ){
if(((i-(N+2))==c(NoSum[-1])[t])) t=t+1
postpred[,i]<-exp(rnorm(iter,K_i[t],1/k_tau[t]^0.5))
}
postpred[,(M+N+2)]<-runif(iter,PO_s,beta)

for (i in 1:length(k_tau)){
j=(M+N+3)+i-1
postpred[,j]<-rgamma(iter,(alpha_ij^2)/alpha_ij_sd^2,(alpha_ij^2)/alpha_ij_sd^2)
}

postpred[,(M+2*N+3)]<-rnorm(iter,r_s,gamma)

for (i in 1:length(r_i)){
j=(M+2*N+4)+i-1
postpred[,j]<-rnorm(iter,r,gamma_i)
}

t=1
for (i in (M+2*N+4):(M+3*N+3) ){
if(((i-(M+2*N+4))==c(NoSum[-1])[t])) t=t+1	
postpred[,i]<-exp(rnorm(iter,r_i[t],gamma_ij))
}

for (i in 1:length(r_tau)){
j=(2*M+3*N+4)+i-1
postpred[,j]<-rgamma(iter,(gamma_ij^2)/gamma_ij_sd^2,(gamma_ij_sd^2)/alpha_ij_sd^2)
}

for (i in 1:length(taui)){
j=(2*M+4*N+4)+i-1
postpred[,j]<-rgamma(iter,(tau^2)/delta^2,tau/delta^2)
}

postpred[,(2*M+5*N+4)]<-rgamma(iter,(tau_s^2)/delta^2,tau_s/delta^2)
postpred
}

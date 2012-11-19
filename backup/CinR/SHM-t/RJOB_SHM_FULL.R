CAPN=6000
load("M_SHM_FULL_27.RData")
library(rjags)
#library(qfa,lib="~/R")
#library(qfaBayes,lib="~/R")
library(qfa,lib="/home/b0919573/R")                                                                                                                        
library(qfaBayes,lib="/home/b0919573/R")  
#################################################
#################################################
priors<-read.table("priors.txt",header=T)

QFA.P[1:18]=priors[1:18,1]


write("
model {
      for (i in 1:N){
      	  for (j in 1:NoORF[i]){
	      	 for (l in 1:NoTime[(NoSum[i]+j)]){
			y[j,l,i] ~ dt(y.hat[j,l,i], exp(nu_l[i]),5)
 			y.hat[j,l,i] <- (K_lm[(NoSum[i]+j)]*P*exp(r_lm[(NoSum[i]+j)]*x[j,l,i]))/(K_lm[(NoSum[i]+j)]+P*(exp(r_lm[(NoSum[i]+j)]*x[j,l,i])-1))
			}
		K_lm[(NoSum[i]+j)]<- exp(K_lm_L[(NoSum[i]+j)])
		K_lm_L[(NoSum[i]+j)] ~ dnorm(K_o_l[i],exp(tau_K_l[i])) I(,0)
		r_lm[(NoSum[i]+j)]<- exp(r_lm_L[(NoSum[i]+j)])
		r_lm_L[(NoSum[i]+j)] ~ dnorm(r_o_l[i],exp(tau_r_l[i]))I(,3.5)
		}
	K_o_l[i] ~ dt( exp(K_p), exp(sigma_K_o),3 ) I(0,)
	r_o_l[i] ~ dt( exp(r_p), exp(sigma_r_o),3 ) I(0,)
	nu_l[i] ~ dnorm(nu_p,  exp(sigma_nu) )

	tau_K_l[i]~dnorm(tau_K_p,sigma_tau_K)I(0,)
	tau_r_l[i]~dnorm(tau_r_p,sigma_tau_r)
	}

K_p ~ dnorm(K_mu,eta_K_p)
r_p ~ dnorm(r_mu,eta_r_p)
nu_p ~ dnorm(nu_mu,eta_nu_p)
P<-exp(P_L)
P_L ~ dnorm(P_mu,eta_P)

sigma_nu~dnorm(eta_nu,psi_nu)
sigma_K_o ~ dnorm(eta_K_o,psi_K_o)
sigma_r_o ~ dnorm(eta_r_o,psi_r_o)

tau_K_p ~ dnorm(tau_K_mu,eta_tau_K_p)
tau_r_p ~ dnorm(tau_r_mu,eta_tau_r_p)
sigma_tau_K ~ dnorm(eta_tau_K,psi_tau_K)
sigma_tau_r ~ dnorm(eta_tau_r,psi_tau_r)

}
","model1.bug")
QFA.I$N=min(QFA.I$N,CAPN)
l<-date()
filename="RCode_SHM"
 jags <- jags.model('model1.bug',
                    data = list('x' = QFA.D$x,
                                'y' = QFA.D$y,
                                'N' = QFA.I$N,
 'NoTime' = QFA.I$NoTime,
 'NoORF' = QFA.I$NoORF,
 'NoSum' = QFA.I$NoSum,
    'tau_K_mu'=QFA.P$sigma_K,                'eta_tau_K_p'=QFA.P$phi_K,
    'tau_r_mu'=QFA.P$sigma_r,                'eta_tau_r_p'=QFA.P$phi_r,
                'eta_tau_K'=QFA.P$phi_K,   'psi_tau_K'=QFA.P$phi_K,
  'eta_tau_r'=QFA.P$phi_r,                 'psi_tau_r'=QFA.P$phi_r,
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
ll<-date()

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
              10000,thin=10)
lll<-date()
save(samp,file=paste(filename,"_F0.R",sep=""))

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
              10000,thin=10)

save(samp,file=paste(filename,"_F1.R",sep=""))
llll<-date()

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
              10000,thin=10)
lllll<-date()
save(samp,file=paste(filename,"_F2.R",sep=""))


stop()

M=4042
pdf(file="testplot2.pdf")
#K_lm[%i]
for (i in 1:M){
j=i
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);
}


#tau_K_l[%i]
j=M+1
for (i in (2*M+3*N+8):(2*M+4*N+7)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);j=j+1
}

#"K_o_l[%i] 
j=M+N+1
for (i in (M+1):(M+N)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);j=j+1
}

#sigma_K_o ");
i=2*M+3*N+5
j=M+2*N+1
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);

#K_p ");
i=M+1+N
j=M+2*N+2
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);

#"P_l ");
i=(M+N+2)
j=M+2*N+3
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);

#r_lm[%i] 
j=M+2*N+4
for (i in (M+2*N+4):(2*M+2*N+3)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);
j=j+1
}

#tau_r_l[%i] ",l);
j=2*M+2*N+4
for (i in (2*M+4*N+8):(2*M+5*N+7)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);
j=j+1
}

#r_o_l[%i] ",l);
j=2*M+3*N+4
for (i in (2*M+2*N+4):(2*M+3*N+3)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);
j=j+1
}

#sigma_r_o ");
i=2*M+3*N+7
j=2*M+4*N+4
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);

#r_p ");
i=2*M+3*N+4
j=2*M+4*N+5
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);

#"nu_l[%i] ",l);
j=2*M+4*N+6
for (i in (M+N+3):(M+2*N+2)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);
j=j+1
}

#sigma_nu ");
i=2*M+3*N+6
j=2*M+5*N+6
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);

#nu_p ");
i=M+2*N+3
j=2*M+5*N+7
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));points(samp[,i],col=2);

dev.off()

i=M+N+3;j=i+N+1;
colnames(samp)[i];names(aa)[j]




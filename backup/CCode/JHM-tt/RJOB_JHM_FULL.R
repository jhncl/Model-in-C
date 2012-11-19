filename="RCode_JHM"
#library(qfaBayes,lib="~/R")
# library(qfa,lib="~/R")
load("M_JHM_FULL_27_27.RData")
library(rjags)
library(qfa,lib="/home/b0919573/R")                                                                                                                        
library(qfaBayes,lib="/home/b0919573/R")  
priors<-read.table("priors.txt",header=T)

QFA.P<-list(
sigma_K=7,		phi_K=0.1,
eta_K_o=8,		psi_K_o=1,
sigma_r=-1,		phi_r=0.1,
eta_r_o=1,		psi_r_o=1,
eta_nu=-1,		psi_nu=1,
K_mu=log(0.2192928),	eta_K_p=1,
r_mu=log(2.5),		eta_r_p=1,
nu_mu=log(31),		eta_nu_p=1,
P_mu=log(0.0002),	eta_P_p=1/0.01,

alpha_mu=0,		eta_alpha=1/(1.5*1.5),
beta_mu=0,		eta_beta=1/(1.5*1.5),
p=0.05,   
eta_gamma=-3.583519,	psi_gamma=1/(4*4),
eta_omega=-3.583519,	psi_omega=1/(4*4)
)
QFA.P[1:18]=priors[1:18,1]
QFA.P[19:27]=priors[20:28,1]


write("
model {
	for (l in 1:N){
		for (c in 1:2){
   		 	for (m in 1:NoORF[l,c]){
	 		        for (n in 1:NoTime[NoSum[l,c]+m,c]){
				y[m,n,l,c] ~ dnorm(y.hat[m,n,l,c],exp(nu_l[l,c]))
				y.hat[m,n,l,c] <- (K_clm[(SHIFT[c]+NoSum[l,c]+m)]*P*exp(r_clm[(SHIFT[c]+NoSum[l,c]+m)]*x[m,n,l,c]))/(K_clm[(SHIFT[c]+NoSum[l,c]+m)]+P*(exp(r_clm[(SHIFT[c]+NoSum[l,c]+m)]*x[m,n,l,c])-1))
				}
			K_clm[(SHIFT[c]+NoSum[l,c]+m)]<-exp(K_clm_L[(SHIFT[c]+NoSum[l,c]+m)])
			K_clm_L[(SHIFT[c]+NoSum[l,c]+m)] ~ dnorm(K_o_l[l]*exp(alpha_c[c]+(delta_l[l,c]*log(gamma_cl[l,c]))),exp(tau_K_cl[l,c]))T(,0)
         		r_clm[(SHIFT[c]+NoSum[l,c]+m)]<-exp(r_clm_L[(SHIFT[c]+NoSum[l,c]+m)])				           
    			r_clm_L[(SHIFT[c]+NoSum[l,c]+m)] ~ dnorm(r_o_l[l]*exp(beta_c[c]+(delta_l[l,c]*log(omega_cl[l,c]))),exp(tau_r_cl[l,c]))T(,3.5)
			}
			tau_K_cl[l,c]~dnorm(tau_K_p[c],exp(sigma_tau_K[c]))T(0,)
			tau_r_cl[l,c]~dnorm(tau_r_p[c],exp(sigma_tau_r[c]))
			nu_l[l,c]~dnorm(nu_p,exp(sigma_nu))
			}
		K_o_l[l]~dt(exp(K_p),exp(sigma_K_o),3)T(0,)
		r_o_l[l]~dt(exp(r_p),exp(sigma_r_o),3)T(0,)
		delta_l[l,1]<-0
                delta_l[l,2]~dbern(p)
		gamma_cl[l,1]<-1
		gamma_cl[l,2]~dt(1,exp(sigma_gamma),3)T(0,)
		omega_cl[l,1]<-1
		omega_cl[l,2]~dt(1,exp(sigma_omega),3)T(0,)
	}
	alpha_c[1]<-0
	alpha_c[2]~dnorm(alpha_mu,eta_alpha)
	beta_c[1]<-0
	beta_c[2]~dnorm(beta_mu,eta_beta)
	
	K_p~dnorm(K_mu,eta_K_p)
	r_p~dnorm(r_mu,eta_r_p)
	nu_p~dnorm(nu_mu,eta_nu_p)
	P <- exp(P_L)
	P_L ~dnorm(P_mu,eta_P)
	sigma_K_o~dnorm(eta_K_o,psi_K_o)
	sigma_r_o~dnorm(eta_r_o,psi_r_o)
	sigma_nu~dnorm(eta_nu,psi_nu)
	sigma_gamma~dnorm(eta_gamma,psi_gamma)
	sigma_omega~dnorm(eta_omega,psi_omega)

tau_K_p[1] ~ dnorm(tau_K_mu,eta_tau_K_p)
tau_K_p[2] ~ dnorm(tau_K_mu,eta_tau_K_p)
tau_r_p[1] ~ dnorm(tau_r_mu,eta_tau_r_p)
tau_r_p[2] ~ dnorm(tau_r_mu,eta_tau_r_p)
sigma_tau_K[1] ~ dnorm(eta_tau_K,psi_tau_K)
sigma_tau_K[2] ~ dnorm(eta_tau_K,psi_tau_K)
sigma_tau_r[1] ~ dnorm(eta_tau_r,psi_tau_r)
sigma_tau_r[2] ~ dnorm(eta_tau_r,psi_tau_r)

}","model1.bug")

library("rjags")
jags <- jags.model('model1.bug',
                   data = list('x' = QFA.D$x,'y' = QFA.D$y,'SHIFT'=QFA.I$SHIFT,
'N' = QFA.I$N,'NoTime' = QFA.I$NoTime,'NoORF' = QFA.I$NoORF,'NoSum' = QFA.I$NoSum, 
'tau_K_mu'=QFA.P$sigma_K,		'eta_tau_K_p'=QFA.P$phi_K,
'tau_r_mu'=QFA.P$sigma_r,		'eta_tau_r_p'=QFA.P$phi_r,
'eta_tau_K'=QFA.P$phi_K,'psi_tau_K'=QFA.P$phi_K,
'eta_tau_r'=QFA.P$phi_K,'psi_tau_r'=QFA.P$phi_K,

'eta_K_o'=QFA.P$eta_K_o,		'psi_K_o'=QFA.P$psi_K_o,
'eta_r_o'=QFA.P$eta_r_o,		'psi_r_o'=QFA.P$psi_r_o,

'eta_nu'=QFA.P$eta_nu,		'psi_nu'=QFA.P$psi_nu,
'K_mu'=QFA.P$K_mu,		'eta_K_p'=QFA.P$eta_K_p,
'r_mu'=QFA.P$r_mu,		'eta_r_p'=QFA.P$eta_r_p,
'nu_mu'=QFA.P$nu_mu,		'eta_nu_p'=QFA.P$eta_nu_p,
'P_mu'=QFA.P$P_mu,		'eta_P'=QFA.P$eta_P,
'alpha_mu'=QFA.P$alpha_mu,		'eta_alpha'=QFA.P$eta_alpha,
'beta_mu'=QFA.P$beta_mu,		'eta_beta'=QFA.P$eta_beta,
'p'=QFA.P$p,   
'eta_gamma'=QFA.P$eta_gamma,	'psi_gamma'=QFA.P$psi_gamma,
'eta_omega'=QFA.P$eta_omega,	'psi_omega'=QFA.P$psi_omega
),
n.chains = 1,n.adapt = 100)
l<-date()
write.table(l,file=paste(filename,"_F",0,"_time.txt",sep=""))
for (i in 1:10){
samp<-coda.samples(jags,
 c(
'K_clm_L',
'tau_K_cl',
'K_o_l',
'sigma_K_o',
'K_p',
'P_L',
'r_clm_L',
'tau_r_cl',
'r_o_l',
'sigma_r_o',
'r_p',
'nu_l',
'sigma_nu',
'nu_p',
'alpha_c',
'beta_c',
'delta_l',
'gamma_cl',
'sigma_gamma',
'omega_cl',
'sigma_omega',
'tau_K_p',
'tau_r_p',
'sigma_tau_K',
'sigma_tau_r'
),
              10000,thin=10)
l<-data()
write.table(l,file=paste(filename,"_F",i,"_time.txt",sep=""))
save(samp,file=paste(filename,"_F",i,".R",sep=""))
}
stop()
save(samp,file=paste(filename,"_F0.R",sep=""))
update(jags,200000)
samp<-coda.samples(jags,
 c(
'K_clm_L',
'tau_K_cl',
'K_o_l',
'sigma_K_o',
'K_p',
'P_L',
'r_clm_L',
'tau_r_cl',
'r_o_l',
'sigma_r_o',
'r_p',
'nu_l',
'sigma_nu',
'nu_p',
'alpha_c',
'beta_c',
'delta_l',
'gamma_cl',
'sigma_gamma',
'omega_cl',
'sigma_omega'
),
              20000,thin=10)
llll<-date()
write.table(file="RJobtime.txt",c(lll,llll))
save(samp,file=paste(filename,"_F1.R",sep=""))
update(jags,200000)
samp<-coda.samples(jags,
 c(
'K_clm_L',
'tau_K_cl',
'K_o_l',
'sigma_K_o',
'K_p',
'P_L',
'r_clm_L',
'tau_r_cl',
'r_o_l',
'sigma_r_o',
'r_p',
'nu_l',
'sigma_nu',
'nu_p',
'alpha_c',
'beta_c',
'delta_l',
'gamma_cl',
'sigma_gamma',
'omega_cl',
'sigma_omega'
),
              2000,thin=10)
lllll<-date()
save(samp,file=paste(filename,"_F2.R",sep=""))

stop()


L=100
M=800
pdf(file="JHMtestplot100.pdf")
#K_clm
t=1
for (i in 1:M)
{
j=i
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
}
dev.off()

#tau_K_cl
j=M+1
for (i in (2*M+9*L+15):(2*M+11*L+14))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#K_o_l
j=M+2*L+1
for (i in (M+1):(M+L))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_K_o
i=2*M+9*L+9
j=M+3*L+1
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#K_p
i=M+L+1
j=M+3*L+2
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#P
i=M+L+2
j=M+3*L+3
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#r_clm
j=M+3*L+4
for (i in (M+8*L+8):(2*M+8*L+7))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#tau_r_cl
j=2*M+3*L+4
for (i in (2*M+11*L+15):(2*M+13*L+14))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#r_o_l
j=2*M+5*L+4
for (i in (2*M+8*L+8):(2*M+9*L+7))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_r_o
i=2*M+9*L+13
j=2*M+6*L+4
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#r_p
i=2*M+9*L+8
j=2*M+6*L+5
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))



#nu_l
j=2*M+6*L+6
for (i in (M+5*L+7):(M+6*L+6))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_nu
i=2*M+9*L+11
j=2*M+7*L+6
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#nu_p
i=M+6*L+7
j=2*M+7*L+7
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))


#alpha_c
i=M+L+4
j=2*M+7*L+8
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#beta_c
i=M+L+6
j=2*M+7*L+9
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#delta_l
j=2*M+7*L+10
for (i in (M+2*L+7):(M+3*L+6))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#gamma_cl
j=2*M+8*L+10
for (i in (M+4*L+7):(M+5*L+6))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_gamma
i=2*M+9*L+10
j=2*M+9*L+10
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))


#omega_cl
j=2*M+9*L+11
for (i in (M+7*L+8):(M+8*L+7))
{
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_omega
i=2*M+9*L+12
j=2*M+10*L+11
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#upsilon_c
i=2*M+13*L+16
j=2*M+10*L+12
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#sigma_upsilon
i=2*M+13*L+16
j=2*M+10*L+13
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

dev.off()

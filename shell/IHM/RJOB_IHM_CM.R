filename="M_IHM_100"
CAPORF=6000
library(qfaBayes)
 library(qfa)
library(rjags)



control<-"dataA2.txt"
query<-"dataB2.txt"


#print(sample)

def="MDRMDP"
upd=10
iter=10 
thin=10

#sample="full"
#work<-paste(upd,iter,thin,sample,sep="_")
#print(work)

#def="CUSTOM"
#REMOVEZEROS<-TRUE

#CUSTOMDEF<-function(x){fun8(x)}
#REMOVEZEROS<-FALSE
#CUT=1000


#################################
print("P")
#################################
p=0.05

#####################################################################
print("Preprocessing")
#####################################################################
a<-read.table(control,header=F)
b<-read.table(query,header=F)
a<-a[,1]
b<-b[,1]
LMN<-read.table("LMNmaxdataA1.txt",header=TRUE)

N<-min(LMN[1,1],CAPORF)
NoORFa<-c(read.table("NoORFdataA1.txt",header=TRUE))$x[1:N]
NoORFb<-c(read.table("NoORFdataB1.txt",header=TRUE))$x[1:N]
M<-sum(NoORFa)

P=exp((a[M+2*N+3]))
K<-exp((a[1:M]))
r<-exp((a[M+2*N+4-1+c(1:M)]))

K[K<=P]=2*P

vecMDR=vecMDP=0
for (i in 1:M){
vecMDR[i]<-r[i]/log(2*max(K[i]-P,0)/max(0,K[i]-2*P)) #MDR
vecMDP[i]<-log(K[i]/P)/log(2) #MDP
}
defa<-vecMDR*vecMDP

M<-sum(NoORFb)

P=exp((b[M+2*N+3]))
K<-exp((b[1:M]))
r<-exp((b[M+2*N+4-1+c(1:M)]))

K[K<=P]=2*P

vecMDR=vecMDP=0
for (i in 1:M){
vecMDR[i]<-r[i]/log(2*max(K[i]-P,0)/max(0,K[i]-2*P)) #MDR
vecMDP[i]<-log(K[i]/P)/log(2) #MDP
}

defb<-vecMDR*vecMDP

NoORF<-cbind(NoORFa,NoORFb)
dimr<-max(NoORF)

NoSuma=NoSumb=0
NoSuma[1]=NoSumb[1]=0
for (i in 2:(1+N)){
NoSuma[i]=NoSuma[i-1]+NoORFa[i-1]
NoSumb[i]=NoSumb[i-1]+NoORFb[i-1]
}

vec=numeric()
for (i in 1:N)
{
vec=c(vec,
c(defa[c(1+NoSuma[i]):NoSuma[i+1]],
rep(NA,dimr-length(defa[c(1+NoSuma[i]):NoSuma[i+1]])),
defb[c(1+NoSumb[i]):NoSumb[i+1]],
rep(NA,dimr-length(defb[c(1+NoSumb[i]):NoSumb[i+1]])))
)
}

y=array(vec,dim=c(dimr,2,N))

################################################
print("Model")
################################################

write("

model {
	for (i in 1:N){
		for (j in 1:2){
			for (k in 1:NoORF[i,j]){
				y[k,j,i]~ dnorm(exp(alpha_c[j]+Z_l[i]+delta_l[i,j]*gamma_cl[i,j]),exp(nu_l[i]+upsilon_c[j]))
			}
		}
		Z_l[i]~dt(Z_p,exp(sigma_Z),4)
		nu_l[i]~dnorm(nu_p,exp(sigma_nu))
		delta_l[i,1]<-0
		delta_l[i,2]~dbern(p)

		gamma_cl[i,1]<-0
		gamma_cl[i,2]~dnorm(0,exp(sigma_gamma))
	}

	alpha_c[1]<-0
	alpha_c[2]~dnorm(alpha_mu,eta_alpha)
	upsilon_c[1]<-0
	upsilon_c[2]~dnorm(upsilon_mu,exp(sigma_upsilon))

	Z_p~dnorm(Z_mu,eta_Z_p)
	nu_p~dnorm(nu_mu,eta_nu_p)
	sigma_Z~dnorm(eta_Z,psi_Z)
	sigma_nu~dnorm(eta_nu,psi_nu)
	sigma_gamma~dnorm(eta_gamma,psi_gamma)
	sigma_upsilon~dnorm(eta_upsilon,psi_upsilon)
}
","model1.bug")

###################################
print("Fit Model")
###################################


priors<-read.table("priors.txt",header=T)
QFA.P=list('Z_mu'=log(50),'eta_Z_p'=1/(6*6),  
'eta_Z'=-6,		'psi_Z'=1/(0.1*0.1),  
'eta_nu'=-2.77,		'psi_nu'=1/(3*3),
'nu_mu'=10.59663,	'eta_nu_p'=1/(5*5), 
'alpha_mu'=0,		'eta_alpha'=1/(1.5*1.5),
'p'=0.05,
'eta_gamma'=-3.583519,          'psi_gamma'=1/(4*4),
'eta_upsilon'=-3.218,		'psi_upsilon'=1,   
'upsilon_mu'=0)
QFA.P[1:16]=priors[1:16,1]

filename="RCode_IHM"

library("rjags")
l<-date()
jags <- jags.model('model1.bug',
data = list('y'=y,
'NoORF'=NoORF,
'N'=N,
'Z_mu'=QFA.P$Z_mu,		'eta_Z_p'=QFA.P$eta_Z_p,  
'eta_Z'=QFA.P$eta_Z,		'psi_Z'=QFA.P$psi_Z,  
'eta_nu'=QFA.P$eta_nu,		'psi_nu'=QFA.P$psi_nu,
'nu_mu'=QFA.P$nu_mu,	'eta_nu_p'=QFA.P$eta_nu_p, 
'alpha_mu'=QFA.P$alpha_mu,		'eta_alpha'=QFA.P$eta_alpha,
'p'=QFA.P$p,
'eta_gamma'=QFA.P$eta_gamma,          'psi_gamma'=QFA.P$psi_gamma,
'eta_upsilon'=QFA.P$eta_upsilon,		'psi_upsilon'=QFA.P$psi_upsilon,   
'upsilon_mu'=QFA.P$upsilon_mu),
n.chains = 1,n.adapt = 100)

###################################
print("Fit Model/UPDATE")
###################################
ll<-date()
update(jags,100000)
lll<-date()
samp<-coda.samples(jags,
          c(
        'Z_l',
        'sigma_Z',
        'Z_p',
        'nu_l',
        'sigma_nu',
        'nu_p',
        'gamma_cl',
        'delta_l',
        'alpha_c',
        'sigma_gamma',
        'upsilon_c',
        'sigma_upsilon'
),
            10000,thin=10)
print("1")
llll<-date()
save(samp,file=paste(filename,"_F0.R",sep=""))
update(jags,100000)
print("2")
print("3")
lllll<-date()
samp<-coda.samples(jags,
          c(
	'Z_l',
	'sigma_Z',     
	'Z_p',  
	'nu_l',          
	'sigma_nu',   
  	'nu_p',  
	'gamma_cl', 
	'delta_l',  
	'alpha_c',
	'sigma_gamma',
	'upsilon_c', 
	'sigma_upsilon'
),
            10000,thin=10)
llllll<-date()
save(samp,file=paste(filename,"_F1.R",sep=""))

stop()
samp<-coda.samples(jags,
          c(
        'Z_l',
        'sigma_Z',
        'Z_p',
        'nu_l',
        'sigma_nu',
        'nu_p',
        'gamma_cl',
        'delta_l',
        'alpha_c',
        'sigma_gamma',
        'upsilon_c',
        'sigma_upsilon'
),
            1000000,thin=10)
lllllll<-date()
save(samp,file=paste(filename,"_F2.R",sep=""))



stop()

plot(samp[,colnames(samp)=="upsilon_c[2]"])

if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-samp}
namesamp<-names(vecsamp)
#write.table(samp,"backup.txt")
#write.table(vecsamp,"backup2.txt")
Z_l<-exp(vecsamp[1:(N)])
sigma_Z<-exp(vecsamp[N+1])
Z<-exp(vecsamp[N+2])
nu_l<-exp(vecsamp[(N+3):(2*N+2)])
sigma_nu<-exp(vecsamp[2*N+3])
nu<-exp(vecsamp[(2*N+4)])
A1<-exp(0)
A2<-exp(vecsamp[2*N+5])
delta<-vecsamp[(2*N+6):(3*N+5)]
gamma<-vecsamp[(3*N+6):(4*N+5)]
sigma_gamma<-exp(vecsamp[(4*N+6)])
upsilon_c <-exp(vecsamp[(4*N+7)])
sigma_upsilon<-exp(vecsamp[(4*N+8)])
if(nrow(samp)>1) {delta_gamma<-colMeans(samp[,(2*N+6):(3*N+5)]*samp[,(3*N+6):(4*N+5)])} else {delta_gamma<-(samp[,(2*N+6):(3*N+5)]*samp[,(3*N+6):(4*N+5)])}
delta_gamma=exp(delta_gamma)

##########################################
print("Ordering")
##########################################
sig<-sum(rep(1,N)[delta>0.5])
order<-order(1-delta)
vecorder<-order(1-delta)[1:sig]

###############
print("Workflow")
###############
pdf(paste("Plots_",work,".pdf",sep=""))
funplot1()
funplot2()
funplot3()
dev.off()
pdf(paste("Plots_Diag_",work,".pdf",sep=""))
funplot4()
funplot5()
dev.off()
pdf(paste("Plots_Indiv_",work,".pdf",sep=""))
funplot6()
dev.off()
funtable()
print(paste("Plots_",work,".pdf",sep=""))
print(paste("Plots_Diag_",work,".pdf",sep=""))
print(paste("Plots_Indiv_",work,".pdf",sep=""))
print(paste("Table_",work,".pdf",sep=""))
#Table<-read.table("Table_terminal1_full.txt",skip=1)
stop()########################################################



stop()

L=100
pdf(file="IHMtestplot100_2.pdf")

#Z_l
for (i in 1:L){
j=i
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])),type="l");plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
}

#sigma_Z
i=6*L+5
j=L+1
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#Z_p
i=L+1
j=L+2
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#nu_l
j=L+3
for (i in (5*L+4):(6*L+3)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_nu
i=6*L+7
j=2*L+3
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#nu_p
i=6*L+4
j=2*L+4
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#alpha_c[2]
i=L+3
j=2*L+5
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#delta_l
j=2*L+6
for (i in (2*L+4):(3*L+3)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#gamma_cl
j=3*L+6
for (i in (4*L+4):(5*L+3)){
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
j=j+1
}

#sigma_gamma
i=6*L+6
j=4*L+6
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
#upsilon_c
i=6*L+10
j=4*L+7
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))

#sigma_upsilon
i=6*L+8
j=4*L+8
plot(density(aa[,j]),main=paste(colnames(samp)[i],t.test((aa[,j]),samp[,i])$p.value));lines(density(samp[,i]),col=2); 
par(mfrow=c(2,1))
plot(c(aa[,j]),main=paste(mean(aa[,j])-mean(samp[,i])));plot(samp[,i],col=2,type="l");
par(mfrow=c(1,1))
dev.off()

i=M+N+3;j=i+N+1;
colnames(samp)[i];names(aa)[j]

###
stop()


load(".RData")
Treat
gene
file=""
samp<-read.table(file,header=T)
if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-samp}
namesamp<-names(vecsamp)
#write.table(samp,"backup.txt")
#write.table(vecsamp,"backup2.txt")
Z_l<-exp(vecsamp[1:(N)])
sigma_Z<-exp(vecsamp[N+1])
Z<-exp(vecsamp[N+2])
nu_l<-exp(vecsamp[(N+3):(2*N+2)])
sigma_nu<-exp(vecsamp[2*N+3])
nu<-exp(vecsamp[(2*N+4)])
A1<-exp(0)
A2<-exp(vecsamp[2*N+5])
delta<-vecsamp[(2*N+6):(3*N+5)]
gamma<-vecsamp[(3*N+6):(4*N+5)]
sigma_gamma<-exp(vecsamp[(4*N+6)])
upsilon_c <-exp(vecsamp[(4*N+7)])
sigma_upsilon<-exp(vecsamp[(4*N+8)])
if(nrow(samp)>1) {delta_gamma<-colMeans(samp[,(2*N+6):(3*N+5)]*samp[,(3*N+6):(4*N+5)])} else {delta_gamma<-(samp[,(2*N+6):(3*N+5)]*samp[,(3*N+6):(4*N+5)])}
delta_gamma=exp(delta_gamma)


pdf(paste("IHM_plot_",file,".pdf",sep=""))

limmin<-0
limmax<-max(A1*Z_l, A2*Z_l*delta_gamma)
i=1:N
plot(1,type="n",main=expression(paste("Treatment, ",Treat,degree,"C"," (delta=Posterior Expectations)")),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control Fitness (=Z_l)",ylab="Query Fitness (=exp(alpha+Z_l+delta_l*gamma_l))",col=8,pch=19,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),col="grey",lwd=2,)
lines(A1*c(-1000,10000),A2*c(-1000,10000),lwd=2,col="grey",lty=2)
lines(c(Z_l[gene=="HIS3"],Z_l[gene=="HIS3"]),c(-1000,10000),lwd=2,col="cadetblue")
lines(c(-1000,10000),c(c(Z_l+delta_gamma)[gene=="HIS3"],c(Z_l+delta_gamma)[gene=="HIS3"]),lwd=2,col="cadetblue")
points(A1*Z_l[i], A2*(Z_l[i]*delta_gamma[i]),col=8,pch=19,cex=0.5)
i=vecorder[delta_gamma[vecorder]>0]
points(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),col=3,pch=19,cex=0.5)
i=vecorder[delta_gamma[vecorder]<=0]  
points(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),col=2,pch=19,cex=0.5)
i=vecorder
text(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),gene[i],pos=4,offset=0.1,cex=0.4)

dev.off()




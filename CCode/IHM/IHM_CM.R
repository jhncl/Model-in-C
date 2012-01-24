filename="M_IHM_100"
CAPORF=4292
library(qfaBayes,lib="~/R")
 library(qfa,lib="~/R")
library(rjags,lib="~/R")



control<-"CCode100Adam.txt"
query<-"CCode100C.txt"


#print(sample)

def="MDRMDP"
upd=10
iter=10 
thin=10

sample="full"
work<-paste(upd,iter,thin,sample,sep="_")
print(work)

#def="CUSTOM"
#REMOVEZEROS<-TRUE

CUSTOMDEF<-function(x){fun8(x)}
REMOVEZEROS<-FALSE
CUT=1000


#################################
print("P")
#################################
p=0.05

#############################################
print("Functions")
##############################################
fun1<-function(x){
as.character(a$Gene[a$ORF%in%x][1])
}
fun2<-function(x){
x[is.na(x)]=-Inf
x[x<0]=NA
x
}
fun3<-function(x){
length(defa[a$ORF==x])
}
fun4<-function(x){
length(defb[b$ORF==x])
}
fun5<-function(x){
c(defa[a$ORF==x],rep(NA,dimr-length(defa[a$ORF==x])),defb[b$ORF==x],rep(NA,dimr-length(defb[b$ORF==x])))
}
fun6<-function(x){
x$r/log(2*(x$K-x$g)/(x$K-2*x$g)) #MDR
}
fun7<-function(x){
log(x$K/x$g)/log(2) #MDP
}
fun8<-function(x){
vecMDRa<-x$r/log(2*(x$K-x$g)/(x$K-2*x$g)) #MDR
vecMDRa<-vecMDPa<-log(x$K/x$g)/log(2) #MDP
as.numeric(lapply(vecMDPa,fun2))*as.numeric(lapply(vecMDRa,fun2))
}
fun9<-function(x){CUSTOMDEF(x)}


########################################
print("Table")
##########################################
funtable<-function(){
alpha1<-A1
alpha2<-A2
nuj1<-nuj[1]
nuj2<-nuj[2]
write.table(rbind(gene,ORF,mu_i,gamma,delta,taui,alpha1,alpha2,mu,tau,nu,nuj1,nuj2)[,order],paste("Table_",work,".pdf",sep=""))
}

###########################################
print("plot fitted with Conditioning on delta=1")
###########################################
funplot1<-function(){
limmin<-min(A1*mu_i, A2*(mu_i+gamma))
limmax<-max(A1*mu_i, A2*(mu_i+gamma))
i=1:N
plot(1,type="n",main=paste("Treatment",treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control (=Alpha1*mu_i)",ylab="Query (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
lines(A1*c(0,1000),A2*c(0,1000),lwd=2)
points(A1*mu_i[i], A2*(mu_i[i]+gamma[i]),main=paste("Treatment",treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single (=Alpha1*mu_i)",ylab="Double (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
i=vecorder[gamma[vecorder]>0]
points(A1*mu_i[i],A2*(mu_i[i]+gamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[gamma[vecorder]<=0]  
points(A1*mu_i[i],A2*(mu_i[i]+gamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(A1*mu_i[i],A2*(mu_i[i]+gamma[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}

###########################################
print("plot fitted NO conditioning")
###########################################
funplot2<-function(){
limmin<-min(A1*mu_i, A2*(mu_i+deltagamma))
limmax<-max(A1*mu_i, A2*(mu_i+deltagamma))
i=1:N
plot(1,type="n",main=paste("Treatment",treat,"Degrees","(delta=Posterior Expectations)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control (=Alpha1*mu_i)",ylab="Query (=Alpha2*(mu_i+delta_i*gamma_i))",col=8,pch=19,cex=0.5)
lines(A1*c(0,1000),A2*c(0,1000),lwd=2)
points(A1*mu_i[i], A2*(mu_i[i]+deltagamma[i]),main=paste("Treatment",treat,"Degrees","(delta=Posterior Expectations)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single (=Alpha1*mu_i)",ylab="Double (=Alpha2*(mu_i+delta_i*gamma_i))",col=8,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]>0]
points(A1*mu_i[i],A2*(mu_i[i]+deltagamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]<=0]  
points(A1*mu_i[i],A2*(mu_i[i]+deltagamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(A1*mu_i[i],A2*(mu_i[i]+deltagamma[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}

##########################################
print("Data plot with highlighted Interactions")
##########################################
funplot3<-function(){
limmin<-min(y,na.rm=TRUE)
limmax<-max(y,na.rm=TRUE)*1.1
plot(1,type="n",main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control",ylab="Query",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),A2*c(-1000,10000),col="cadetblue",lwd=2,)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=2)
abline(lm(colMeans(y[,2,],na.rm=TRUE)~0+colMeans(y[,1,],na.rm=TRUE)),col="grey",lty=3)
lines(c(mean(defa[gene=="HIS3"]),mean(defa[gene=="HIS3"])),c(-1000,1000),lwd=2)
lines(c(-1000,1000),c(mean(defb[gene=="HIS3"]),mean(defb[gene=="HIS3"])),lwd=2)
points(colMeans(y[,1,],na.rm=TRUE),colMeans(y[,2,],na.rm=TRUE),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,col=8,cex=0.5)
i=vecorder[deltagamma[vecorder]>0]
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]<=0]  
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),gene[i],pos=4,offset=0.1,cex=0.4)
legend(1,limmax, c("Model Fit (y=Alpha2*x)","1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("cadetblue","grey","grey","black"), lty=c(1,2,3,1))

plot(1,type="n",main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control",ylab="Query",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),A2*c(-1000,10000),col="cadetblue",lwd=2,)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=2)
abline(lm(colMeans(y[,2,],na.rm=TRUE)~0+colMeans(y[,1,],na.rm=TRUE)),col="grey",lty=3)
lines(c(mean(defa[gene=="HIS3"]),mean(defa[gene=="HIS3"])),c(-1000,1000),lwd=2)
lines(c(-1000,1000),c(mean(defb[gene=="HIS3"]),mean(defb[gene=="HIS3"])),lwd=2)
points(colMeans(y[,1,],na.rm=TRUE),colMeans(y[,2,],na.rm=TRUE),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,col=8,cex=0.5)
i=vecorder[deltagamma[vecorder]>0]
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]<=0]  
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=1:N
text(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),gene[i],pos=4,offset=0.1,cex=0.4)
legend(1,limmax, c("Model Fit (y=Alpha2*x)","1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("cadetblue","grey","grey","black"), lty=c(1,2,3,1))
}

###########################################
print("Diagnostics trace acf density")
###########################################
funplot4<-function(){
par(mfrow=c(4,1))
den<-matrix(N,2000,8)
den[,1]<-rgamma(2000,(alpha_a^2)/(alpha_b^2),alpha_a/(alpha_b^2))
den[,2]<-rnorm(2000,0,gam_b^(-0.5))
den[,3]<-rnorm(2000,mu_a,mu_b^(-0.5))
den[,4]<-rnorm(2000,mu_a,mu_b^(-0.5))
den[,5]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
den[,6]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
den[,7]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
den[,8]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
namesampden<-unique(substring(namesamp,1,4))
names<-namesampden[-2]
for (i in 1:ncol(den))
{
plot(density(den[,i]),paste(names[i],"Prior Density"))
}
par(mfrow=c(4,3))
for (i in c(2,(3*N+3):c(7*N+5)))
{
plot(as.numeric(samp[,i]),main=paste("Trace Top",namesamp[i],"(Treatment",treat,"Degrees)"),type="l")
plot(density(as.numeric(samp[,i])),main=paste("Density",namesamp[i]))
acf(as.numeric(samp[,i]),main=paste("ACF",namesamp[i]))
}
par(mfrow=c(4,2))
for (i in (N+3):(2*N+2)){
plot(as.numeric(samp[,i]),main=paste("Trace Top",namesamp[i],"(Treatment",treat,"Degrees)"),type="l")
plot(density(as.numeric(samp[,i])),main=paste("Density",namesamp[i]))
}
par(mfrow=c(1,1))
}

#####################################################################
print("Number of Repeats for (Single+Double)Versus Order of Sig Inter")
#####################################################################
funplot5<-function(){
plot(rowSums(NoORF)[order],main="No. of Rep for (Single+Double) Versus Order of Sig Inter",ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",pch=19,cex=0.5)
points(rowSums(NoORF)[order][1:sig],ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",col=2,pch=19,cex=0.5)
plot(rowSums(NoORF)[order],main="No. of Rep for (Single+Double) Versus Order of Sig Inter",ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",pch=19,cex=0.5)
points(rowSums(NoORF)[order][1:sig],ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",col=2,pch=19,cex=0.5)
i=order
text(rowSums(NoORF)[order],gene[i],pos=4,offset=0.1,cex=0.4)
}

##########################################
print("Individual plots")
##########################################
funplot6<-function(){
limmin<-min(y,na.rm=TRUE)
limmax<-max(y,na.rm=TRUE)
vec<-order(1-delta)[1:sig];col=1
for (i in vec)
{
plot(y[,1,i],y[,2,i],main=paste(gene[i],"Interacting","Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,cex=0.5,col=i)
lines(c(-1000,10000),A2*c(-1000,10000),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",cex=0.5,col=i,lwd=2)
}
vec<-order(1-delta)[N:(N-sig)];col=8
for (i in vec)
{
plot(y[,1,i],y[,2,i],main=paste(gene[i],"Non Interacting","Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,cex=0.5,col=i)
lines(c(-1000,10000),A2*c(-1000,10000),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",lwd=2,cex=0.5,col=i)
}
}

#####################################################################
print("Preprocessing")
#####################################################################
a<-read.table(control,header=TRUE)
b<-read.table(query,header=TRUE)
LMN<-read.table("LMNmaxdata.txt",header=TRUE)

N<-min(LMN[1,1],CAPORF)
NoORFa<-c(read.table("NoORFdataA1.txt",header=TRUE))$x[1:N]
NoORFb<-c(read.table("NoORFdataB1.txt",header=TRUE))$x[1:N]
M<-sum(NoORFa)

P=exp(mean(a[,M+2*N+3]))
K<-exp(colMeans(a[1:M]))
r<-exp(colMeans(a[M+2*N+4-1+c(1:M)]))

vecMDR=vecMDP=0
for (i in 1:M){
vecMDR[i]<-r[i]/log(2*max(K[i]-P,0)/max(0,K[i]-2*P)) #MDR
vecMDP[i]<-log(K[i]/P)/log(2) #MDP
}
defa<-vecMDR*vecMDP
M<-sum(NoORFb)

P=exp(mean(b[,M+2*N+3]))
K<-exp(colMeans(b[1:M]))
r<-exp(colMeans(b[M+2*N+4-1+c(1:M)]))

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

vec=0
for (i in 1:N)
{
vec=c(vec,
c(defa[c(1+NoSuma[i]):NoSuma[i+1]],
rep(NA,dimr-length(defa[c(1+NoSuma[i]):NoSuma[i+1]])),
defb[c(1+NoSumb[i]):NoSumb[i+1]],
rep(NA,dimr-length(defb[c(1+NoSumb[i]):NoSumb[i+1]])))
)
}

vec<-vec[-1]


y=array(vec,dim=c(dimr,2,N))


#######################################
#print("Data Storage")
#######################################
#write.table(y,"y_i.txt")
#write.table(NoORF,"NoORF_i.txt")
#NoORF<-read.table("NoORF_i.txt")
#NoORF<-as.matrix(NoORF)
#N=nrow(NoORF)
#y<-read.table("y_i.txt")
#y<-array(c(as.matrix(y)),dim=c(max(NoORF),2,N))

#########################################
print("Statistics")
#########################################
#var(na.omit(c(y[,1,])))
#var(na.omit(c(y[,2,])))

###################################
print("TESTING")
###################################
#yy<-y
#y<-y[,,1:sample]
#N=sample


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
		Z_l[i]~dnorm(Z_p,exp(sigma_Z))
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
sigma_Z<-min(10,sigma_Z_UT)
	sigma_Z_UT~dnorm(eta_Z,psi_Z)
	sigma_nu~dnorm(eta_nu,psi_nu)
	sigma_gamma~dnorm(eta_gamma,psi_gamma)
	sigma_upsilon~dnorm(eta_upsilon,psi_upsilon)
}
","model1.bug")

###################################
print("Fit Model")
###################################
library("rjags")
jags <- jags.model('model1.bug',
data = list('y'=y,
'NoORF'=NoORF,
'N'=N,
'Z_mu'=log(50),		'eta_Z_p'=1/(6*6),  
'eta_Z'=-2.772589,	'psi_Z'=1/(7*7),  
'eta_nu'=-2.77,		'psi_nu'=1/(3*3),
'nu_mu'=10.59663,	'eta_nu_p'=1/(5*5), 
'alpha_mu'=0,		'eta_alpha'=1/(1.5*1.5),
'p'=0.05,
'eta_gamma'=-3.583519,          'psi_gamma'=1/(4*4),
'eta_upsilon'=-3.218,		'psi_upsilon'=1,   
'upsilon_mu'=0
),
n.chains = 1,n.adapt = 100)

###################################
print("Fit Model/UPDATE")
###################################
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
            1000,thin=1)
print("1")
save(samp,file=paste(filename,"_F0.R",sep=""))

print("2")

update(jags, 10000)
print("3")
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

save(samp,file=paste(filename,"_F1.R",sep=""))


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

save(samp,file=paste(filename,"_F2.R",sep=""))



stop()

plot(samp[,colnames(samp)=="upsilon_c[2]"])

if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-samp}
namesamp<-names(vecsamp)
#write.table(samp,"backup.txt")
#write.table(vecsamp,"backup2.txt")
A1<-vecsamp[1]
A2<-vecsamp[2]
delta<-vecsamp[(N+3):(2*N+2)]
gamma<-vecsamp[(3*N+3):(4*N+2)]
if(nrow(samp)>1) {deltagamma<-colMeans(samp[,(N+3):(2*N+2)]*samp[,(3*N+3):(4*N+2)])} else {deltagamma<-(samp[,(N+3):(2*N+2)]*samp[,(3*N+3):(4*N+2)])}
mu_i<-vecsamp[(4*N+4):(5*N+3)]
mu<-vecsamp[4*N+3]
nu<-vecsamp[(5*N+4)]
nuj<-vecsamp[(5*N+5):(5*N+6)]
tau<-vecsamp[(5*N+7)]
taui<-vecsamp[(5*N+8):(6*N+7)]

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




library(qfaBayes,lib="~/R")
 library(qfa,lib="~/R")
 Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")

DescripControl<-"ExptDescriptionCDC13.txt"
a<-rod.read(files=Control[1],inoctimes=DescripControl)
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))
Treat<-27
MPlate<-15#as.character(unique(a$MasterPlate.Number))
a<-funcREMOVE(a,Screen,Treat,MPlate)

Row<-a$Row
Col<-a$Col
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
IDstrip=sapply(IDstrip,unique)
IDstrip=lapply(IDstrip,i=6,funcStrip)
IDstrip=unlist(IDstrip)
IDstrip=na.omit(IDstrip)
a<-a[a$ID%in%IDstrip,]
#########

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]


Control<-c("cdc13-1_rad9D_SDLv2_Rpt1.txt","cdc13-1_rad9D_SDLv2_Rpt2.txt","cdc13-1_rad9D_SDLv2_Rpt3.txt","cdc13-1_rad9D_SDLv2_Rpt4.txt")
DescripControl<-"ExptDescriptionCDC13RAD9.txt"
b<-rod.read(files=Control[1],inoctimes=DescripControl)

qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))
Treat<-27
MPlate<-15#as.character(unique(b$MasterPlate.Number))
b<-funcREMOVE(b,Screen,Treat,MPlate)


Row<-b$Row
Col<-b$Col
for (i in 1:nrow(b)){
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

b$ID<-paste(b$Barcode,b$MasterPlate.Number,Row,Col,sep="")
ORFuni_b<-unique(b$ORF)
####
sum(rep(1,length(ORFuni))[ORFuni==ORFuni_b])/length(ORFuni)
####

########
funcIDlist<-function(x){
b$ID[b$ORF==x]
}
funcStrip<-function(x,i){x[1:i]}
IDstrip=sapply(ORFuni,funcIDlist)
IDstrip=sapply(IDstrip,unique)
IDstrip=lapply(IDstrip,i=6,funcStrip)
IDstrip=unlist(IDstrip)
IDstrip=na.omit(IDstrip)
b<-b[b$ID%in%IDstrip,]
#########

b<-b[order(b$ORF,b$ID,b$Expt.Time), ]

#a<-funcIDORDER(a)
IDuni<-unique(a$ID)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))

#b<-funcIDORDER(b)
IDuni<-unique(b$ID)
###



#gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);M<-length(IDuni)#?
NoORF_b<-unlist(lapply(ORFuni,funcNoORF,data=b))#no of repeats each orf
NoTime_b<-c(0,unlist(lapply(IDuni,funcNoTime,data=b)))# 0+ no of time each repeat
NoSum_b<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_b)))


dimr<-max(NoORF_a,NoORF_b);dimc<-max(NoTime_a,NoTime_b)

y<-funcXY_J(a$Growth,b$Growth,M,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)
x<-funcXY_J(a$Expt.Time,b$Expt.Time,M,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)

QFA.I<-list("NoORF"=cbind(NoORF_a,NoORF_b),"NoTime"=cbind(NoTime_a,NoTime_b)[-1,],"NoSum"=cbind(NoSum_a,NoSum_b),"N"=N,"M"=M,"gene"=gene,SHIFT=c(0,max(NoSum_a,NoSum_b))
)
Scaling=TRUE######
if (Scaling==TRUE){y<-funcSCALING(rbind(a,b),y)}
QFA.D<-list(x=x,y=y)


Priors<-list(
sigma_K=13,		phi_K=3,
eta_K_o=13,		psi_K_o=3,
sigma_r=-1,		phi_r=3,
eta_r_o=13,		psi_r_o=3,
eta_nu=15,		psi_nu=1,
K_mu=log(0.2192928),	eta_K_p=1,
r_mu=log(2.5),		eta_r_p=1,
nu_mu=log(31),		eta_nu_p=1,
P_mu=log(0.0002),	eta_P_p=1/0.00001,

alpha_mu=log(1),	eta_alpha=1/3^2,
beta_mu=log(1),		eta_beta=1/4^2,
p=0.05,   
eta_gamma=1,	psi_gamma=1,
eta_omega=1,	psi_omega=1,
eta_upsilon=1,	psi_upsilon=1,	    
upsilon_mu=1,	eta_upsilon_p=1)

write("
model {
	for (i in 1:N){
		for (c in 1:2){
   		 	for (j in 1:NoORF[i,c]){
	 		        for (l in 1:NoTime[NoSum[i,c]+j,c]){
				y[j,l,i,c] ~ dnorm(y.hat[j,l,i,c],exp(upsilon_c[c]+nu_l[l]))
				y.hat[j,l,i,c] <- (K_clm[(SHIFT[c]+NoSum[i,c]+j)]*P*exp(r_clm[(SHIFT[c]+NoSum[i,c]+j)]*x[j,l,i,c]))/(K_clm[(SHIFT[c]+NoSum[i,c]+j)]+P*(exp(r_clm[(SHIFT[c]+NoSum[i,c]+j)]*x[j,l,i,c])-1))
				}
			K_clm[(SHIFT[c]+NoSum[i,c]+j)]<-exp(K_clm_L[(SHIFT[c]+NoSum[i,c]+j)])
			K_clm_L[(SHIFT[c]+NoSum[i,c]+j)] ~ dnorm(alpha_c[c]+(K_o_l[i]+delta_l[i,c]*gamma_cl[i,c]),exp(tau_K_cl[i+(c-1)*N]))
         		r_clm[(SHIFT[c]+NoSum[i,c]+j)]<-exp(min(3.5,r_clm_L[(SHIFT[c]+NoSum[i,c]+j)]))				           
    			r_clm_L[(SHIFT[c]+NoSum[i,c]+j)] ~ dnorm(beta_c[c]+(r_o_l[i]+delta_l[i,c]*omega_cl[i,c]),exp(tau_r_cl[i+(c-1)*N]))
			}
			tau_K_cl[i+(c-1)*N]~dnorm(sigma_K,phi_K)
			tau_r_cl[i+(c-1)*N]~dnorm(sigma_r,phi_r)
		}

		K_o_l[i] ~ dnorm(K_p,exp(sigma_K_o))
		r_o_l[i] ~ dnorm(r_p,exp(sigma_r_o))
		delta_l[i,1]<-0
                delta_l[i,2]~dbern(p)
		nu_l[i]~ dnorm(nu_p,exp(sigma_nu))
		gamma_cl[i,1]<-0
		gamma_cl[i,2]~dnorm(0,exp(sigma_gamma))
		omega_cl[i,1]<-0
		omega_cl[i,2]~dnorm(0,exp(sigma_omega))
	}
	alpha_c[1]<-1
	alpha_c[2]~dnorm(alpha_mu,eta_alpha)
	beta_c[1]<-1
	beta_c[2]~dnorm(beta_mu,eta_beta)
	upsilon_c[1]~dnorm(upsilon_p,exp(sigma_upsilon))
	upsilon_c[2]~dnorm(upsilon_p,exp(sigma_upsilon))

	
	K_p ~ dnorm(K_mu,eta_K_p)
	r_p ~ dnorm(r_mu,eta_r_p)
	nu_p~dnorm(nu_mu,eta_nu_p)
	upsilon_p~dnorm(upsilon_mu,eta_upsilon_p)
	P <- exp(P_L)
	P_L ~ dnorm(P_mu,eta_P)
	sigma_K_o~dnorm(eta_K_o,psi_K_o)
	sigma_r_o~dnorm(eta_r_o,psi_r_o)
	sigma_nu~dnorm(eta_nu,psi_nu)
	sigma_upsilon~dnorm(eta_upsilon,psi_upsilon)
	sigma_gamma~dnorm(eta_gamma,psi_gamma)
	sigma_omega~dnorm(eta_omega,psi_omega)

}","model1.bug")





 QFA.P<-Priors
library("rjags")
jags <- jags.model('model1.bug',
                   data = list('x' = QFA.D$x,'y' = QFA.D$y,'SHIFT'=QFA.I$SHIFT,
'N' = QFA.I$N,'NoTime' = QFA.I$NoTime,'NoORF' = QFA.I$NoORF,'NoSum' = QFA.I$NoSum, 
'sigma_K'=QFA.P$sigma_K,		'phi_K'=QFA.P$phi_K,
'eta_K_o'=QFA.P$eta_K_o,		'psi_K_o'=QFA.P$psi_K_o,
'sigma_r'=QFA.P$sigma_r,		'phi_r'=QFA.P$phi_r,
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
'eta_omega'=QFA.P$eta_omega,	'psi_omega'=QFA.P$psi_omega,
'eta_upsilon'=QFA.P$eta_upsilon,	'psi_upsilon'=QFA.P$psi_upsilon,	    
'upsilon_mu'=QFA.P$upsilon_mu,	'eta_upsilon_p'=QFA.P$eta_upsilon_p
),n.chains = 1,n.adapt = 100)
date()
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
'upsilon_c',
'sigma_upsilon',
'upsilon_p'
),
              1000,thin=1)

save(samp,file="MTEST2_0.R")

samp<-coda.samples(jags,
 c(
'lnK_ij',
'lnr_ij',
'K_ij',
'r_ij',
'tau_K_l',                                                                      'tau_r_l',                                           
'gamma',  
'omega',  
'delta',
'K_o_l',
'r_o_l',
'alpha',
'beta',
'nu_c',
'nu',
'lnP',
'K_p',
'r_p',
'upsilon_l',
'upsilon_p',
'sigma_Ko',
'sigma_ro',
'sigma_nu',
'sigma_upsilon',
'sigma_gamma',
'sigma_omega'
),
              100000,thin=100)

save(samp,file="MTEST2_1.R")

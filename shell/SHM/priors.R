
#load("M_SHM_FULL.RData")

#P#
t=vec_P=1
for (i in 1:QFA.I$N){
    for (j in 1:QFA.I$NoORF[i]){
    	vec_P[t]<-QFA.D$y[j,1,i]
		t=t+1
		}
}
P=median(vec_P)

#K_lm#
t=l=1
for (i in 1:QFA.I$N){
    for (j in 1:QFA.I$NoORF[i]){
    	l[t]<-QFA.D$y[j,QFA.I$NoTime[ QFA.I$NoSum[i]+j ],i]
				      t=t+1
				      }
}
K_lm=l


#K_l#
t=l=1
for (i in 1:QFA.I$N){
l[t]=median(K_lm[QFA.I$NoSum[i]:QFA.I$NoSum[i+1]])
t=t+1
}
K_l=l



#K#
K=median(K_l)

#r#
t=l=1
for (i in 1:QFA.I$N){
    for (j in 1:QFA.I$NoORF[i]){
    	val<-order(na.omit(QFA.D$y[j,-1,i]-QFA.D$y[j,-QFA.I$NoTime[QFA.I$NoSum[i]+j],i]),decreasing=T)[1]
		valt<-(QFA.D$x[j,val+1,i]+QFA.D$x[j,val,i])/2
		valt2<-(QFA.D$y[j,val+1,i]+QFA.D$y[j,val,i])/2
			l[t]<-log(valt2*(K_lm[t]-P)/(P*(K_lm[t]-valt2)))/valt
				t=t+1
				}
}
r_lm=l
r_lm[is.nan(r_lm)]=2.5
r_lm[r_lm<0.01][K_lm[r_lm<0.01]>(P)]=0.01
r_lm[r_lm<0.01]=0.01

#r_l#
t=l=1
for (i in 1:QFA.I$N){
l[t]=median(r_lm[QFA.I$NoSum[i]:QFA.I$NoSum[i+1]])
t=t+1
}
r_l=l

#r#
r=median(r_l)
if(r==0){r=mean(r_l);if(r==0){r=0.000000000000000001}}
#sigma_K_o#
l=K_l
l[l<=0]=min(l[!l<=0])
sigma_K_o=1/var(log(l))

#tau_K_l#
t=l=ll=1
for (i in 1:QFA.I$N){
ll=K_lm[QFA.I$NoSum[i]:QFA.I$NoSum[i+1]]
ll[ll==0]=min(ll[!(ll==0)])
l[t]<-1/var(log(ll))
t=t+1
}
tau_K_l=l
tau_K_l[tau_K_l==Inf]=max(tau_K_l[!tau_K_l==Inf])
tau_K_l[is.na(tau_K_l)]=max(tau_K_l[!is.na(tau_K_l)])

#sigma_r_o#
l=r_l
l[l<=0]=min(l[!l<=0])
sigma_r_o=1/var(log(l))

#tau_r_l#
t=l=1
for (i in 1:QFA.I$N){
ll=r_lm[QFA.I$NoSum[i]:QFA.I$NoSum[i+1]]
ll[ll==0]=min(ll[!(ll==0)])
l[t]<-1/var(log(ll))
t=t+1
}
tau_r_l=l
tau_r_l[tau_r_l==Inf]=max(tau_r_l[!tau_r_l==Inf])
tau_r_l[is.na(tau_r_l)]=max(tau_r_l[!is.na(tau_r_l)])

#nu_l#
t=l=ll=1
for (i in 1:QFA.I$N){
    for (j in 1:QFA.I$NoORF[i]){
    	ll[t]<-var(QFA.D$y[j,c(1,2),i])
		t=t+1
		}
}
t=llll=1
for (i in 1:QFA.I$N){
llll[t]=1/median(ll[QFA.I$NoSum[i]:QFA.I$NoSum[i+1]])
t=t+1
}
nu_l=llll
nu_l[nu_l==Inf]=max(nu_l[!nu_l==Inf])
#nu_l#
nu_mu=(nu_l)

l=nu_l
l[l==0]=min(l[!l==0])
sigma_nu=1/var(log(l))

###
###
 

QFA.P<-list(
 sigma_K=log(median(tau_K_l)),  phi_K=1/(4*abs(log(mean(tau_K_l))))^2,#gd
 eta_K_o=log(sigma_K_o), psi_K_o=1/(12*abs(log(mean(sigma_K_o))))^2,#gd
 sigma_r=log(mean(tau_r_l)),  phi_r=1/(4*abs(log(mean(tau_r_l))))^2,#gd
 eta_r_o=log(sigma_r_o), psi_r_o=1/(2*abs(log(mean(sigma_r_o))))^2,#gd

 eta_nu=log(sigma_nu),    psi_nu=1/(10*abs(log(median(sigma_nu))))^2,#gd

 K_mu=log(K),		              eta_K_p=1/(3*abs(log(median(K))))^2,
 r_mu=log(r),               eta_r_p=1/(3*abs(log(median(r))))^2,
 nu_mu=log(median(nu_mu)),          eta_nu_p=1/(2*abs(log(median(nu_mu))))^2,
 P_mu=log(P),        eta_P=1/(2*(abs(log(median(P))))^2),
N=c(QFA.I$N)
)
write.table(as.matrix(QFA.P),file=paste("priors_",Treat,".txt",sep=""))
write.table(cbind(tau_K_l,tau_r_l,nu_l,QFA.I$NoORF),file=paste("para_est_l_",Treat,".txt",sep="")  )
write.table(cbind(K_lm,r_lm,vec_P),file=paste("para_est_lm_",Treat,".txt",sep="")  )

#library(gplots,lib="/home/b0919573/R/x86_64-pc-linux-gnu-library/2.10")

pdf(paste("priors_SHM_",Treat,".pdf",sep="")) 
par(mfcol=c(2,1))
t=1
a=QFA.P$sigma_K
bb=1/(QFA.P$phi_K)^0.5
l<-rnorm(2000,a,bb)
plot(density(l),main=names(QFA.P)[t])#log
lines(density(log(tau_K_l)),col=2) #log
textplot(data.frame(quantile(l),quantile(log(tau_K_l))))

t=t+2
a=QFA.P$eta_K_o
bb=1/(QFA.P$psi_K_o)^0.5
l<-rnorm(2000,a,bb)
plot(density(l),main=names(QFA.P)[t])
lines(c(QFA.P$eta_K_o,QFA.P$eta_K_o),c(0,10000),col=2) #log
lines(c(QFA.P$sigma_K,QFA.P$sigma_K),c(0,10000),col=3) #log
textplot(data.frame(quantile(l),QFA.P$eta_K_o,QFA.P$sigma_K))

t=t+2
a=QFA.P$sigma_r
bb=1/(QFA.P$phi_r)^0.5
l<-rnorm(2000,a,bb)
plot(density(l),main=names(QFA.P)[t])#log
lines(density(log(tau_r_l)),col=2) #log
textplot(data.frame(quantile(l),quantile(log(tau_r_l))))


t=t+2
a=QFA.P$eta_r_o
bb=1/(QFA.P$psi_r_o)^0.5
l<-rnorm(2000,a,bb)
plot(density(l),main=names(QFA.P)[t])
lines(c(QFA.P$eta_r_o,QFA.P$eta_r_o),c(0,10000),col=2) #log
lines(c(QFA.P$sigma_r,QFA.P$sigma_r),c(0,10000),col=3) #log
textplot(data.frame(quantile(l),QFA.P$eta_r_o,QFA.P$sigma_r))

t=t+2
a=QFA.P$eta_nu
bb=1/(QFA.P$psi_nu)^0.5
l<-rnorm(2000,a,bb)
plot(density(l),main=names(QFA.P)[t])
lines(c(QFA.P$eta_nu,QFA.P$eta_nu),c(0,10000),col=2) #log
textplot(data.frame(quantile(l),QFA.P$eta_nu))


t=t+2
a=QFA.P$K_mu
bb=1/(QFA.P$eta_K_p)^0.5
plot(density((rnorm(2000,a,bb))),main=names(QFA.P)[t])  #log
lines(density(log(K_l)),col=2) #log
textplot(data.frame((quantile((rnorm(2000,a,bb)))),quantile(log(K_l))))#log

t=t+2
a=QFA.P$r_mu
bb=1/(QFA.P$eta_r_p)^0.5
plot(density((rnorm(2000,a,bb))),main=names(QFA.P)[t])
lines(density(log(r_l)),col=2) #log
textplot(data.frame((quantile((rnorm(2000,a,bb)))),quantile(log(r_l))))#log

t=t+2
a=QFA.P$nu_mu
bb=1/(QFA.P$eta_nu_p)^0.5
l<-rnorm(2000,a,bb)
plot(density(l),main=names(QFA.P)[t])
lines(density(log(nu_l)),col=2) #log
textplot(data.frame(quantile(l),QFA.P$eta_nu))

t=t+2
a=QFA.P$P_mu
bb=1/(QFA.P$eta_P)^0.5
plot(density((rnorm(2000,a,bb))),main=names(QFA.P)[t])
lines(density(log(vec_P)),col=2) #log
textplot(data.frame(quantile(rnorm(2000,a,bb)),quantile(log(vec_P))))#log
t=t+2
dev.off()

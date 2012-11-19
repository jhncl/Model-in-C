load("LARGE.RData")

library(qfa,lib="~/R")
 library(qfaBayes,lib="~/R")
Scaling=TRUE
IDuni<-unique(a$ID)
ORFuni<-unique(a$ORF)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of 
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))
dimr<-max(NoORF_a);dimc<-max(NoTime_a)
y<-funcXY(a$Growth,M,N,NoTime_a,NoSum_a,dimr,dimc)
x<-funcXY(a$Expt.Time,M,N,NoTime_a,NoSum_a,dimr,dimc)

QFA.I<-list("NoORF"=c(NoORF_a),"NoTime"=c(NoTime_a)[-1],"NoSum"=c(NoSum_a),"N"=N,"M"=M,"gene"=gene)
if (Scaling==TRUE){y<-funcSCALING(a,y)}
QFA.D<-list(y=y,x=x,ORFuni=ORFuni)


source("C_CM2.Priors")
 QFA.P<-Priors
 source("C_CM2")

QFA.I$N=1000
 jags <- jags.model('model1.bug',
                    data = list('x' = QFA.D$x,
                                'y' = QFA.D$y,
                                'N' = QFA.I$N,
 'NoTime' = QFA.I$NoTime,
 'NoORF' = QFA.I$NoORF,
 'NoSum' = QFA.I$NoSum,
 'K_s' = QFA.P$K_s,
 'r_s' = QFA.P$r_s,
 'PO_s' = QFA.P$PO_s,
 'beta' = QFA.P$beta,
 'tau_s' = QFA.P$tau_s,
 'delta'=QFA.P$delta,
 'alpha' = QFA.P$alpha,
 'alpha_i' = QFA.P$alpha_i,
 'alpha_i_tau'=QFA.P$alpha_i_tau,
 'alpha_ij' = QFA.P$alpha_ij,
 'alpha_ij_tau'=QFA.P$alpha_ij_tau,
 'gamma' = QFA.P$gamma,
 'gamma_i' = QFA.P$gamma_i,
 'gamma_i_tau'=QFA.P$gamma_i_tau,
 'gamma_ij' = QFA.P$gamma_ij,
 'gamma_ij_tau'=QFA.P$gamma_ij_tau,
 'delta_mu'=QFA.P$delta_mu,
 'delta_sd'=QFA.P$delta_sd
 ),
                    n.chains = 1,
                    n.adapt = 100
,inits=list('K'=log(0.1),'r'=log(1.7),'PO_L'=log(0.0037),'K_i_tau_L'=log(5),'r_i_tau_L'=log(0.02))
)

 
 samp<-coda.samples(jags,
 c('K_ij',
 'r_ij',
 'K_i',
 'r_i',
 'K',
 'PO',
 'r',
 'tau',
 'tau_i',
 'K_ij_tau',
 'r_ij_tau',
 'K_i_tau',
 'r_i_tau',
 'delta_tau'),
              1000,thin=1)
save(samp,file=paste("G_CM3_",0,".R",sep=""))
 for (i in 1:50)
 {
 samp<-coda.samples(jags,
 c('K_ij',
 'r_ij',
 'K_i',
 'r_i',
 'K',
 'PO',
 'r',
 'tau',
 'tau_i',
 'K_ij_tau',
 'r_ij_tau',
 'K_i_tau',
 'r_i_tau',
 'delta_tau'),
              100000,thin=100)
 save(samp,file=paste("G_CM3_",i,".R",sep=""))
}

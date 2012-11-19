#### Hierachical Logistic Curve Model ####
qfa.Hierachical<-function(experiment,Scaling,iter,upd,thin,inits,PlotOutput=TRUE,work,CustomModel=FALSE){
a<-experiment
a<-funcIDORDER(a)
IDuni<-unique(a$ID)
ORFuni<-unique(a$ORF)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))
dimr<-max(NoORF_a);dimc<-max(NoTime_a)
y<-funcXY(a$Growth,M,N,NoTime_a,NoSum_a,dimr,dimc)
x<-funcXY(a$Expt.Time,M,N,NoTime_a,NoSum_a,dimr,dimc)
QFA.I<-list("NoORF"=c(NoORF_a),"NoTime"=c(NoTime_a)[-1],"NoSum"=c(NoSum_a),"N"=N,"M"=M,"gene"=gene)
if (Scaling==TRUE){y<-funcSCALING(a,y)}
QFA.D<-list(y=y,x=x,ORFuni=ORFuni)
if (!(CustomModel==FALSE)){source(CustomModel)} else {funcMODELHierarchical()}
QFA.P<-funcPRIORS(CustomModel)

#if (LinearGaussian==FALSE){funcPriorPlot(QFA.P)} else {funcPriorPlot_LG(QFA.P)} 
#print("Please check Priors.pdf. Are you happy with your choice of #priors? (y/n) ")
#ans<-readline()
#if (!(ans=="y")) stop() else print("Fitting Model")
 
samp<-funcFITandUPDATE(QFA.I,QFA.D,QFA.P,inits,iter,upd,thin)
QFA.O<-funcPosterior(samp,N,M,iter,thin,upd)
QFA<-c(QFA.O,QFA.I,QFA.D,QFA.P)
if(PlotOutput==TRUE){qfaplots.H(QFA,work)}
return(QFA)
}

### Hierachical Logistic Curve Model Plots to Pdf###
qfaplots.H<-function(QFA,work,LinearGaussian=FALSE){

samp<-QFA$samp
iter<-QFA$iter
thin<-QFA$thin

y<-QFA$y
x<-QFA$x

N<-QFA$N
M<-QFA$M
NoSum<-QFA$NoSum
NoORF<-QFA$NoORF
NoTime<-QFA$NoTime
gene<-QFA$gene

K_s<-QFA$K_s
r_s<-QFA$r_s
PO_s<-QFA$PO_s
beta<-QFA$beta
tau_s<-QFA$tau_s
delta<-QFA$delta

alpha<-QFA$alpha
gamma<-QFA$gamma

alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_i_tau<-QFA$alpha_i_tau
gamma_i_tau<-QFA$gamma_i_tau

alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij
alpha_ij_tau<-QFA$alpha_ij_tau
gamma_ij_tau<-QFA$gamma_ij_tau

delta_mu<-QFA$delta_mu
delta_sd<-QFA$delta_sd

namesamp<-QFA$namesamp
K<-QFA$K
K_i<-QFA$K_i
K_ij<-QFA$K_ij
PO<-QFA$PO
r<-QFA$r
r_i<-QFA$r_i
r_ij<-QFA$r_ij
taui<-QFA$taui
tau<-QFA$tau
K_i_tau<-QFA$K_i_tau
r_i_tau<-QFA$r_i_tau
K_ij_tau<-QFA$K_ij_tau
r_ij_tau<-QFA$r_ij_tau

delta_tau<-QFA$delta_tau

if(LinearGaussian==TRUE){
K<-exp(QFA$K)
K_i<-exp(QFA$K_i)
K_ij<-QFA$K_ij
r<-exp(QFA$r)
r_i<-exp(QFA$r_i)
tau<-exp(QFA$tau)
}

################################################
print("priors")
################################################

if (LinearGaussian==FALSE){funcPriorPlot(QFA)} else {funcPriorPlot_LG(QFA)} 




################################################
print("Plots")
################################################

ylimmin<-0
ylimmax<-max(na.omit(as.numeric(y)))
xlimmin<-0
xlimmax<-max(na.omit(as.numeric(x)))

pdf(paste("Plots_M_",work,".pdf",sep=""))
################################################
print("Master Curve")
################################################
plot(x,y,main="Master Curve",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
curve((K*PO*exp(r*x))/(K+PO*(exp(r*x)-1)), xlimmin, xlimmax,add=TRUE,col=1)
################################################
print("ORF Curves")
################################################
plot(-1,-1,main="ORF Curves",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
for (i in 1:N){
points(x[,,i],y[,,i],main="ORF Curves",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax),col=i)
}
for (i in 1:N)
{
curve((K_i[i]*PO*exp(r_i[i]*x))/(K_i[i]+PO*(exp(r_i[i]*x)-1)), 0, 8,add=TRUE,col=i) 
}
################################################
print("Repeat Curves")
################################################
plot(x,y,main="Repeat Curves", xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
for (i in 1:M)
{
curve((K_ij[i]*PO*exp(r_ij[i]*x))/(K_ij[i]+PO*(exp(r_ij[i]*x)-1)), 0, 8,add=TRUE,col=i) 
}


################################################
print("Model Variation tau") #post pred master
################################################
plot(x,y,main="Curve variation tau_m", xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
if (LinearGaussian==FALSE){MCurveVar<-funcMCurveVar(QFA)} else {MCurveVar<-funcMCurveVar_LG(QFA)} 
KK=K
rr=r
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)),xlimmin, xlimmax,add=TRUE) 
lines(MCurveVar$Time,MCurveVar$MQQU,col=2)
lines(MCurveVar$Time,MCurveVar$MQQD,col=2)

dev.off()

pdf(paste("Plots_M_indiv_",work,".pdf",sep=""))
###########################################
print("plots for individual Logistic curve fits")#postpredORF andorfrep
###########################################

if (LinearGaussian==FALSE){ICurveVar<-funcICurveVar(QFA)} else {ICurveVar<-funcICurveVar_LG(QFA)} 

dev.off()

pdf(paste("Plots_Tau_",work,".pdf",sep=""))

plot(x,y,main="Master Curve",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
curve((K*PO*exp(r*x))/(K+PO*(exp(r*x)-1)), xlimmin, xlimmax,add=TRUE)
curve(+2/tau^0.5+(K*PO*exp(r*x))/(K+PO*(exp(r*x)-1)), xlimmin, xlimmax,add=TRUE,col=2)
curve(-2/tau^0.5+(K*PO*exp(r*x))/(K+PO*(exp(r*x)-1)), xlimmin, xlimmax,add=TRUE,col=2)

for (i in 1:N)
{
for (j in 1:NoORF[i])
{
plot(-1,-1,main=paste(gene[i],"Repeat Curves"),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[j,,i],y[j,,i])
KK=K_ij[(j+NoSum[i])]
rr=r_ij[(j+NoSum[i])]
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=1) 
curve(2/taui[i]^0.5+(KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=1+j) 
curve(-2/taui[i]^0.5+(KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=1+j) 
}
}

dev.off()


pdf(paste("Plots_M_diag_",work,".pdf",sep=""))
###########################################
print("Prior density")
###########################################
#par(mfrow=c(4,2))
sampsize<-round(iter/thin)
if (LinearGaussian==FALSE){den<-funcDen(sampsize,QFA)}else den<-funcDen_LG(sampsize,QFA)
namesampden<-unique(substring(namesamp,1,4))
#for (i in 1:ncol(den))
#{
#plot(density(den[,i]),paste(namesampden[i],"Prior Density"))
#}

###########################################
print("Diagnostics trace acf density")
###########################################
if (LinearGaussian==FALSE){postpred<-funcPostPred(sampsize,QFA)} else postpred<-funcPostPred_LG(sampsize,QFA)
################################################

namesamp_EDIT=namesamp=colnames(samp)
namesamp_EDIT[(M+N+3):(M+2*N+2)]<-"K_izj"
namesamp_EDIT[(2*M+3*N+7):(2*M+4*N+6)]<-"r_izj"
namesampden<-unique(substring(namesamp_EDIT,1,4))
#################################################
par(mfrow=c(4,5))
for (i in 1:ncol(samp))
{
post<-density(as.numeric(samp[,i]))
pred<-density(postpred[,i])

plot(as.numeric(samp[,i]),main=paste(namesamp[i],"Trace Top"),type="l")
acf(as.numeric(samp[,i]),main=paste(namesamp[i],"ACF"))

plot(post,main=paste(namesamp[i],"Posterior Density (Black)"),xlim=c(min(post$x),max(post$x)),
ylim=c(min(post$y,post$y),max(post$y,post$y)))
lines(pred,col=2)

plot(post,main=paste(namesamp[i],"Posterior Predictive Density (Red)"),xlim=c(min(pred$x),max(pred$x)),
ylim=c(min(pred$y,pred$y),max(pred$y,pred$y)))
lines(pred,col=2)

t<-(1:ncol(den))[namesampden==substring(namesamp_EDIT[i],1,4)]
pri<-density(den[,t])
plot(pri,main=paste(namesamp[i],"Top level Prior Density"),xlim=c(min(pri$x,pri$x),max(pri$x,pri$x)),
ylim=c(min(pri$y,pri$y),max(pri$y,pri$y)),col=3)

lines(post)
}

dev.off()
}

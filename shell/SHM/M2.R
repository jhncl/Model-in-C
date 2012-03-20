a=aa=aaa=list(0)
t=1
folder_list=dir()[file.info(dir())$isdir]
folder_list=folder_list[!is.na(folder_list)]
for (i in 1:length(folder_list)){
setwd(paste("~/QFADatasets/SHM/",folder_list[i],sep=""))
folder_list_inner=dir()[file.info(dir())$isdir]
folder_list_inner=folder_list_inner[!is.na(folder_list_inner)]
for (j in 1:length(folder_list_inner)){
setwd(paste("~/QFADatasets/SHM/",folder_list[i],"/",folder_list_inner[j],sep=""))
list=list.files()
list=list[grepl("^priors.*\\.txt", list)]
list=list[!c(list=="priors.txt")]
print(list)
list2=list.files()
list2=list2[grepl("^para_est_l_.*\\.txt", list2)]
print(list2)
list3=list.files()
list3=list3[grepl("^para_est_lm_.*\\.txt", list3)]
print(list3)
a[[t]]=read.table(list,header=T)
aa[[t]]=read.table(list2,header=T)
aaa[[t]]=read.table(list3,header=T)
names(a[[t]])=paste(folder_list[i],"_",folder_list_inner[j],sep="")
names(aa[[t]])=paste(folder_list[i],"_",folder_list_inner[j],sep="")
names(aaa[[t]])=paste(folder_list[i],"_",folder_list_inner[j],sep="")
t=t+1
}
setwd("~/QFADatasets/SHM/")
}



LEN=t-1

vec=matrix(nrow=nrow(a[[1]]),ncol=LEN)
for (j in 1:nrow(a[[1]])){
for (i in 1:LEN){
vec[j,i]=a[[i]][j,1]
if (a[[i]][19,1]<1000){vec[j,i]=NA}
}
}


vec2=matrix(nrow=nrow(a[[1]]),ncol=6)
for (i in 1:nrow(a[[1]]))
{
vec2[i,]=c(quantile( vec[i,][!is.na(vec[i,])]),mean(vec[i,][!is.na(vec[i,])]))
}

vec3=numeric(0)
for (i in 2*c(1:9)-1){
vec3[i]=vec2[i,3]
vec3[i+1]=1/( (max(vec2[i,5]-vec2[i,3],vec2[i,3]-vec2[i,1])*1.1)/2 )^2
}

l=numeric(0)
for (i in 1:LEN){
l=c(l,log(aa[[i]][,1]/2))}### t dist only df =4 ###########################

ll=numeric(0)
for (i in 1:LEN){
ll=c(ll,log(aa[[i]][,2]))}
lll=numeric(0)
for (i in 1:LEN){
lll=c(lll,log(aa[[i]][,3]))}

vec4=a[[1]]
names(vec4)="MASTER"
vec4[,1]=c(vec3,"N")
vec4[1,1]=median(l)#sigma_K 
vec4[2,1]=1/( (max(abs(median(l)-max(l)),median(l)-min(l))*1.1)/2 )^2
vec4[5,1]=median(ll)#sigma_r 
vec4[6,1]=1/( (max(abs(median(ll)-max(ll)),median(ll)-min(ll))*1.1)/2 )^2
vec4[15,1]=median(lll)#nu_mu
vec4[16,1]=1/( (max(abs(median(lll)-max(lll)),median(lll)-min(lll))*1.1)/2 )^2

write.table(as.numeric(vec4[,1]),file="priors.txt")
vec5=rbind(vec4,
alpha_mu=0,             eta_alpha=1/(2)^2,
beta_mu=0,              eta_beta=1/(2)^2,
p=0.05,
eta_gamma=vec4[3,1],     psi_gamma=vec4[4,1],
eta_omega=vec4[7,1],    psi_omega=vec4[8,1],
eta_upsilon=log(1/(2*2)),     psi_upsilon=1/(2)^2,
upsilon_mu=0)
write.table(as.numeric(vec5[,1]),file="priors_JHM.txt")

pdf(file="priors.pdf")
for (i in seq(1,nrow(a[[1]]),2) ){
plot(density(vec[i,][!is.na(vec[i,])]),main=c(vec2[i,3],vec2[i,6]))
}
dev.off()


stop()

K=exp(as.numeric(vec4[11,1]))
r=exp(as.numeric(vec4[13,1]))
P=exp(as.numeric(vec4[17,1]))
Z_mu=log((log((K/P))/log(2))*(r/log(2*((K-P)/(K-2*P)))))
K_vec=exp(range(vec[11,][!is.na(vec[11,])]))
r_vec=exp(range(vec[13,][!is.na(vec[13,])]))
l=log((log((K_vec/P))/log(2))*(r_vec/log(2*((K_vec-P)/(K_vec-2*P)))))
eta_Z_p=1/(max(abs(Z_mu-l[1]),abs(Z_mu-l[2])) /2 )^2
eta_alpha=1/((l[2]/(l[1]))/2)^2
p=0.05

K=r=P=numeric(0)
for (i in 1:LEN){
r=c(r,(aaa[[i]][,2]))
P=c(P,(aaa[[i]][,3]))
K=c(K,aaa[[i]][,1])
}
P[P<=0]=min(P[!(P<=0)])
K[K<2*P]=2*P[K<2*P]
r[r<=0]=0
r[r>8]=8
vec=(log((K/P))/log(2))*(r/log(2*((K-P)/(K-2*P))))
NoORF=numeric(0)
for (i in 1:LEN){
NoORF=c(NoORF,aa[[i]][,4])}

l=numeric(0)
t=1
for (i in 1:length(NoORF)){
l[i]=1/var(vec[t:c(t+NoORF[i]-1)])
t=t+NoORF[i]
}
l=log(l)
nu_mu=median(l)
l=l[!is.na(l)]
eta_nu_p=1/((max(abs(median(l)-max(l)),median(l)-min(l))*1.1)/2 )^2
#
ll=numeric(0)
for (i in 1:LEN){
ll[i]=nrow(aaa[[i]])
}

Z=ZZ=numeric(0)
t=1
for (i in 1:length(NoORF)){
Z[i]=median(vec[t:c(t+NoORF[i]-1)])
t=t+NoORF[i]
}

t=1
for (i in 1:LEN){
ZZ[i]=1/var(Z[t:c(t-1+ll[i])])
}
ZZ=log(ZZ)
eta_Z=median(ZZ)
psi_Z=1/((max(abs(median(ZZ)-max(ZZ)),median(ZZ)-min(ZZ))*1.1)/2 )^2


t=1
lll=numeric(0)
for (i in 1:LEN){
lll[i]=1/var(is.finite(l[t:c(t-1+a[[i]][19,1])]))
t=t+a[[i]][19,1];i=i+1
}
lll=log(lll)
lll[is.na(lll)]=mean(lll[!is.na(lll)])
eta_nu=median(lll)
psi_nu=1/((max(abs(median(lll)-max(lll)),median(lll)-min(lll))*1.1)/2 )^2


#1/( (max(abs(median(l)-max(l)),median(l)-min(l))*1.1)/2 )^2

QFA.P=list(
'Z_mu'=Z_mu,	'eta_Z_p'=eta_Z_p,  
'eta_Z'=eta_Z,		    'psi_Z'=psi_Z,  
'eta_nu'=eta_nu,		    'psi_nu'=psi_nu,
'nu_mu'=nu_mu,			    'eta_nu_p'=eta_nu_p, 
'alpha_mu'=0,			    'eta_alpha'=eta_alpha,
'p'=p,
'eta_gamma'=eta_Z,		'psi_gamma'=psi_Z,   
'eta_upsilon'=log(1/(2*2)),     'psi_upsilon'=1/(2)^2,
'upsilon_mu'=0
)
write.table(as.numeric(QFA.P),file="priors_IHM.txt")

K=exp(as.numeric(vec4[11,1]))
r=exp(as.numeric(vec4[13,1]))
P=exp(as.numeric(vec4[17,1]))
Z_mu=((log((K/P))/log(2))*(r/log(2*((K-P)/(K-2*P)))))
Z_mu=log(1.01^Z_mu-1)
K_vec=exp(range(vec[11,][!is.na(vec[11,])]))
r_vec=exp(range(vec[13,][!is.na(vec[13,])]))
l=((log((K_vec/P))/log(2))*(r_vec/log(2*((K_vec-P)/(K_vec-2*P)))))
l=log(1.01^l-1)
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
vec=log(1.01^vec-1)
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

nu_mu=median(l[is.finite(l)])

eta_nu_p=1/((max(abs(median(l[is.finite(l)])-max(l[is.finite(l)])),median(l[is.finite(l)])-min(l[is.finite(l)]))*1.1)/2 )^2
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
for (i in 1:LEN){##
ZZ[i]=1/var(Z[t:c(t-1+ll[i])][is.finite(Z[t:c(t-1+ll[i])])])###
}###

ZZg=log(ZZ)
ZZg[is.na(ZZg)]=mean(ZZg[!is.na(ZZg)])
eta_gamma=(median(ZZg))
psi_gamma=1/((max(abs((median(ZZg))-(max(ZZg))),abs((median(ZZg))-(min(ZZg))))*1.1)/2 )^2

ZZ=log(ZZ/2)###### t dist only df =4 ###################################
ZZ[is.na(ZZ)]=mean(ZZ[!is.na(ZZ)])
eta_Z=(median(ZZ))
psi_Z=1/((max(abs((median(ZZ))-(max(ZZ))),abs((median(ZZ))-(min(ZZ))))*1.1)/2 )^2

t=1
lll=numeric(0)
for (i in 1:LEN){
lll[i]=1/var(  l[t:c(t-1+a[[i]][19,1])][ is.finite(l[t:c(t-1+a[[i]][19,1])])])
t=t+a[[i]][19,1];i=i+1
}
lll=log(lll)
lll[is.na(lll)]=mean(lll[!is.na(lll)])
eta_nu=median(lll)
psi_nu=1/((max(abs(median(lll)-max(lll)),median(lll)-min(lll))*1.1)/2 )^2


#1/( (max(abs(median(l)-max(l)),median(l)-min(l))*1.1)/2 )^2

QFA.P=list(
'Z_mu'=Z_mu,    'eta_Z_p'=eta_Z_p,
'eta_Z'=eta_Z,              'psi_Z'=psi_Z,
'eta_nu'=eta_nu,                    'psi_nu'=psi_nu,
'nu_mu'=nu_mu,                      'eta_nu_p'=eta_nu_p,
'alpha_mu'=0,                       'eta_alpha'=eta_alpha,
'p'=p,
'eta_gamma'=eta_gamma,              'psi_gamma'=psi_gamma,
'eta_upsilon'=log(1/(2*2)),     'psi_upsilon'=1/(2)^2,
'upsilon_mu'=0
)
write.table(as.numeric(QFA.P),file="priors_IHM.txt")

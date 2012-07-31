# Hierarchical.R
# Single Model
CustomModel="C_CM2"

Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")
DescripControl<-"ExptDescriptionCDC13.txt"
source("iter.R")
#upd=8
#iter=2
#thin=2
work=paste(upd,iter,thin,sep="_")
data("AdamFull8Rep")
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))
Treat<-27
MPlate<-unique(a$MasterPlate.Number)
a<-funcREMOVE(a,Screen,Treat,MPlate)

inits=list('K'=log(0.089),'r'=log(1.82),'PO_L'=log(0.00224),'K_i_tau_L'=1.25,'r_i_tau_L'=3.429)


a<-a[!(a$Col==24),]
a<-a[!(a$Col==1),]
a<-a[!(a$Row==1),]
a<-a[!(a$Row==16),]

vec<-a[a$ORF=="YOR202W",][1:399,]#peter 392 adam 399
a<-a[!(a$ORF=="YOR202W"),]
a<-rbind(a,vec)

a<-funcIDORDER(a)
IDuni<-unique(a$ID)
t=0
for (i in 1:length(IDuni)){
vec<-a[a$ID==IDuni[i],]$Growth
LENV<-length(vec)
for (j in 2:LENV){
if (vec[j]<vec[j-1]){t=1}
}
if (t==1){a<-a[!a$ID==IDuni[i],]}
t=0
print(i*100/N)
}

ControlFit<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=FALSE,work="ModelHExample",CustomModel=CustomModel,inits=inits)
save(ControlFit,file=paste(CustomModel,work,"R",sep="."))

### Plots ###
qfaplots.H(ControlFit,CustomModel,LinearGaussian=TRUE)

# eof



# Hierarchical.R
# Single Model
CustomModel="A_CustomModel7_2"
work="_L"
Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")
DescripControl<-"ExptDescriptionCDC13.txt"
upd=1000000
iter=1000000
thin=100
data("AdamFull")

qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)

ControlFit<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=FALSE,work="ModelHExample",CustomModel=CustomModel)
save(ControlFit,file=paste(CustomModel,work,"R",sep="."))

### Plots ###
qfaplots.H(ControlFit,CustomModel,LinearGaussian=TRUE)

# eof



# Joint.R
# Joint Model
upd=1000000
iter=1000000
thin=100
work="Joint_L"
PlotOutput=FALSE
data("Adam_cdc-1_SDLV2_REP1")
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))[1]
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)


data("cdc13-1_rad9D_SDLv2_Rpt1")

qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))[1]
b<-funcREMOVE(b,Screen,Treat,MPlate)
CustomModel="A_JCMNNormals2"
JointFit<-qfa.Joint(a,b,Scaling=TRUE,iter,upd,thin,PlotOutput=FALSE,work,CustomModel=FALSE)
save(ControlFit,file=paste(CustomModel,work,"R",sep="."))

qfaplots.J(JointFit,work)

# eof



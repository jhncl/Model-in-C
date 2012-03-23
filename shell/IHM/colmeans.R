
aa<-read.table("CC_SHM_URA3_Raw_27_5000_400000_10000_10.R",header=T)                                     
dataA2=colMeans(aa)                                                        
write.table(dataA2,file="A2_CC_SHM_URA3_Raw_27_5000_400000_10000_10.R",row.names=F,col.names=F)                   

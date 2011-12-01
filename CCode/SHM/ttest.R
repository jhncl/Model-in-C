vec=vec2=0
for( i in 1:ncol(samp))
{
vec[i]=t.test((aa[,i]),samp[,i])$p.value
if(vec[i]<0.05){vec2[i]=
t.test(exp(aa[,i]),samp[,i])$p.value
}
}
vec2[is.na(vec2)]=0

vec=vec2=0
for( i in 1:ncol(samp))
{
vec[i]=var.test((aa[,i]),samp[,i])$p.value
if(vec[i]<0.05){vec2[i]=
var.test(exp(aa[,i]),samp[,i])$p.value
}
}
vec2[is.na(vec2)]=0


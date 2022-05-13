survival.mle <-
function(x,fit,param.vec,censoring,censor.type="right"){

a<-sum(log(dgl(x[which(censoring==1)],fit,param=param.vec)))

if(censor.type=="right"){
b<-sum(log(1-pgl(x[which(censoring==0)],fit,param=param.vec)))}

if(censor.type=="left"){
b<-sum(log(pgl(x[which(censoring==0)],fit,param=param.vec)))}

return(a+b)
}

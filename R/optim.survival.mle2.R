optim.survival.mle2 <-
function(value, censoring,x, y, param,intercept,censor.type){
value<-c(intercept,value)
len<-length(value)
resid <- y - x %*% value[-c((len - 3):len)]
gld.val<-value[c((len - 3):len)]
return(survival.mle(resid,gld.val,param,censoring,censor.type)*-1)

}

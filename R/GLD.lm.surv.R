GLD.lm.surv <-
function (formula, censoring, data, param, maxit = 20000, fun, method = "Nelder-Mead", 
    diagnostics = TRUE, range = c(0.01, 0.99), init = NULL, alpha = 0.05,censor.type="right",adj.int=TRUE,GLD.adj=FALSE,adj.censor=TRUE,keep.uncen=TRUE) 
{

# Reorder based on censoring:

data<-data[order(censoring),]
censoring<-censoring[order(censoring)]

    
    init.mod<-lm(formula, data=data[which(censoring==1),])
    # The only difference from fun.lm.surv.ext1 is that the earlier work use GLD regression in the first step.
   
    len.init<-length(init.mod$coeff)+4
    id.censoring<-which(censoring==1)
    
   
    init.full<-model.frame(formula, data)
    y <- init.full[,1]
    x <- model.matrix(formula,data)
    resid1<-y-x %*% init.mod$coeff
    resid1.fit<-fun(resid1)

    value <- c(init.mod$coeff, resid1.fit[1:4])

    result <- optim(value, optim.survival.mle1, x = x, y = y, censoring=censoring, param = param, censor.type=censor.type,control = list(maxit = 20000), method = method)

    AIC.full<-2*length(result$par)+2*result$value
    BIC.full<-length(result$par)*log(nrow(data))+2*result$value
    
    full.result <- result$par[-c((length(result$par) - 3):length(result$par))]

    if(adj.int==TRUE){

    if(adj.censor==TRUE){
    adj <- mean(y[id.censoring] - x[id.censoring,,drop=F] %*% full.result)}

    # Could also use GLD.lm, the following will lead to lower AIC, but not a desired fit

    if(adj.censor==FALSE){
    adj <- lm(y ~ x %*% full.result)$coeff[1]
    }

    if (is.element("(Intercept)", dimnames(x)[[2]]) == TRUE) {
        full.result[1] <- full.result[1] + adj
    }
    if (is.element("(Intercept)", dimnames(x)[[2]]) == FALSE) {
        warning("Bias adjustment is provided separately from estimated \n parameters, adjustment is necessary to ensure residuals sum to zero")
    }

    resid2.fit<-fun(y - x %*% full.result)
    value<-c(full.result,resid2.fit)  

    result<-optim(value[-1], optim.survival.mle2, x = x, y = y, intercept=value[1],censoring=censoring, censor.type=censor.type,param = param, control = list(maxit = 20000), method = method)

    AIC.full<-2*(length(result$par)-1)+2*result$value
    BIC.full<-(length(result$par)-1)*log(nrow(data))+2*result$value

    full.result <- c(value[1],result$par[-c((length(result$par) - 3):length(result$par))])

    }
    
    if (result$convergence == 0) {
        converge.report <- "converged"
    }
    if (result$convergence != 0) {
        converge.report <- "not converged"
    }
    message1 <- paste("This analysis was carried out using", 
        toupper(param), "GLD")
    message2 <- paste("The error distribution was estimated using", 
        "Maximum likelihood estimation")
    message3 <- paste("The optimisation procedure used was", 
        as.character(substitute(method)), "and it has", converge.report)
    messages <- rbind(message1, message2, message3)
    dimnames(messages)[[1]] <- NULL

    fitted <- x %*% full.result

    if(adj.int==TRUE){
    if (is.element("(Intercept)", dimnames(x)[[2]]) == TRUE) {
        resid <- y - fitted
    }
    if (is.element("(Intercept)", dimnames(x)[[2]]) == FALSE) {
        resid <- y - (fitted + adj)
    }}

   if(adj.int==FALSE){
        adj=0
        resid <- y - fitted
    
    }


    if (GLD.adj==TRUE){
    full.result <- c(full.result, fun.mean.convert(c(result$par[c((length(result$par) - 
        3):length(result$par))]), param))}

    if (GLD.adj==FALSE){
    full.result <- c(full.result, c(result$par[c((length(result$par) - 
        3):length(result$par))]))}

    names(full.result)[1:(length(full.result) - 4)] <- dimnames(x)[[2]]
    names(full.result)[(length(full.result) - 3):length(full.result)] <- paste("L", 
        1:4, sep = "")
 

 gld.values <- full.result[(length(full.result) - 3):length(full.result)]
 ddst.pval<-ddst.uniform.test(pgl(resid[id.censoring],gld.values,param=param),compute.p=T)$p.value
censor.gld.values<-gld.values 

if(ddst.pval<0.05){
full.result[(length(full.result) - 3):length(full.result)]<-fun(resid[id.censoring])
if(keep.uncen==TRUE){
censor.gld.values<-full.result[(length(full.result) - 3):length(full.result)] }
}

    if (diagnostics) {
        gld.values <- full.result[(length(full.result) - 3):length(full.result)]
        par(mfrow = c(1, 1))
	# Range is only active in summaryGraphics.gld.surv.lm for the full version
        qqgld.default(resid[id.censoring], gld.values, param = param)
        legend("bottomright", c(paste(toupper(param), "GLD"), 
            paste("(", paste(signif(gld.values, 3), collapse = ","), 
                ")", sep = "")), bty = "n")
        pval <- ks.gof(resid[id.censoring], "pgl", lambda1 = gld.values, param = param)$p.value
        diag.real <- fun.diag.ks.g(data = resid[id.censoring], result = gld.values, 
            param = param, alpha = alpha)/1000 * 100
   ddst.pval<-ddst.uniform.test(pgl(resid[id.censoring],gld.values,param=param),compute.p=T)$p.value
        legend("topleft", paste("KS test p-value=", format.pval(pval), 
            "\n", "Data driven smooth test p-value=",format.pval(ddst.pval), "\n",
            "Resample KS test > ", alpha, " is", diag.real, "%"), bty = "n")
    }
    print(messages)
    print(full.result)

    return(list(Message = messages, `Bias Correction` = adj, 
        `Estimated parameters` = full.result, Fitted = fitted, 
        Residual = resid, formula = formula, param = param, y = y, 
        x = x, fun = fun,censoring=censoring,AIC.full=AIC.full,BIC.full=BIC.full,censor.gld.values=censor.gld.values))
}

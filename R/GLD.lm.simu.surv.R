GLD.lm.simu.surv <-
function (formula, censoring,data, init.coeff, init.resid, param, maxit = 20000, fun,
    method = "Nelder-Mead", init = NULL,fit=NULL,censor.type="right",adj.int = TRUE, GLD.adj = FALSE, adj.censor = TRUE)
{
    
    value <- c(init.coeff, init.resid)
    fit<<-fit
    censoring<<-censoring
    y <- model.frame(formula, data = data)[,1]
    x <- model.matrix(formula, data = data)

    max.y<-max(na.omit(y))
    censor.id<-which(censoring==0)
    id.censoring<-which(censoring==1)

    if(tolower(censor.type)=="right"){   
    y[censor.id]<-log(sapply(1:length(censor.id),function(i) runif(1,0,exp(y[censor.id][i]))))}

    if(tolower(censor.type)=="left"){
    y[censor.id]<-log(sapply(1:length(censor.id),function(i) runif(1,exp(y[censor.id][i]),exp(max.y))))}

    check <- optim.survival.mle1(value, censoring, x, y, param,censor.type)

    if ((is.na(check) | is.inf(check))) {
        repeat {
            init.mod <- lm(formula, data = data)
            y <- init.mod$model[, 1]
            max.y<-max(na.omit(y))

            if(tolower(censor.type)=="right"){   
    y[censor.id]<-log(sapply(1:length(censor.id),function(i) runif(1,0,exp(y[censor.id][i]))))}

    if(tolower(censor.type)=="left"){
    y[censor.id]<-log(sapply(1:length(censor.id),function(i) runif(1,exp(y[censor.id][i]),exp(max.y))))}

            check <- optim.survival.mle1(value, censoring, x, y, param,censor.type)

            if (is.null(init)) {
                init.1 <- init.mod$coeff
                init.2 <- fun(init.mod$resid)
            }
            else if (!is.null(init) & length(init.mod$coeff) == 
                length(init)) {
                init.1 <- init
                init.2 <- fun(y - x %*% init)
            }
            else if (!is.null(init) & length(init.mod$coeff) != 
                length(init)) {
                init.len <- length(init)
                init.1 <- init[1:(init.len - 4)]
                init.2 <- init[(init.len - 3):init.len]
            }
            value <- c(init.1, init.2)
   
            check <- optim.survival.mle1(value, censoring, x, y, param,censor.type)

            if (is.na(check) == FALSE & is.inf(check) == FALSE) {
                break
            }
        }
    }

    rm(fit)
    result <- optim(value, optim.survival.mle1, x = x, y = y, censoring=censoring, param = param, censor.type=censor.type,control = list(maxit = 20000), method = method)

    full.result <- result$par[-c((length(result$par) - 3):length(result$par))]
    adj<-0

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

    full.result <- c(value[1],result$par[-c((length(result$par) - 3):length(result$par))])

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

    if (result$convergence == 0) {
        converge.report <- "converged"
    }
    if (result$convergence != 0) {
        converge.report <- "not converged"
    }
    message1 <- paste(toupper(param), "GLD")
    message2 <- converge.report
    r <- cbind(adj, t(full.result), message1, message2)
    rm(censoring)
    return(r)
}

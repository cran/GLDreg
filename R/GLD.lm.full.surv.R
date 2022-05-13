GLD.lm.full.surv <-
function (formula, censoring, data, param, maxit = 20000, fun, method = "Nelder-Mead", 
    range = c(0.01, 0.99), n.simu = 1000, summary.plot = FALSE, 
    init = NULL, alpha = 0.05,censor.type="right",adj.int=FALSE,GLD.adj=FALSE,adj.censor=TRUE,keep.uncen=TRUE)  
{
    fit <- GLD.lm.surv(formula, censoring, data, param, maxit, fun, method, diagnostics = FALSE, 
        range, init,alpha,censor.type,adj.int,GLD.adj,adj.censor,keep.uncen)
    formula.info <- unlist(strsplit(deparse(formula), " ~ "))
    if (length(formula.info) == 2 & formula.info[[2]] == ".") {
        data <- na.omit(data)
        fit.simu.f <- fun.simu.gld.lm.alt.surv(n.simu, c(fit$y[which(censoring==1)] + rgl(nrow(fit$Fitted[which(censoring==1),,drop=F]), 
            fit$"Estimated parameters"[match(c("L1", "L2", "L3", 
                "L4"), names(fit$"Estimated parameters"))], param = fit$param),fit$y[which(censoring==0)] + rgl(nrow(fit$Fitted[which(censoring==0),,drop=F]), 
            fit$censor.gld.values, param = fit$param))~., fit, censoring, data = data[, -c(match(formula.info[1], dimnames(data)[[2]]))], 
            param = fit$param, fun = fun, init = init,fit=fit,censor.type=censor.type,adj.int = adj.int, GLD.adj = GLD.adj, adj.censor = adj.censor)
    }
    if (length(formula.info) != 2 | formula.info[[2]] != ".") {
        fit.simu.f <- fun.simu.gld.lm.alt.surv(n.simu, update(fit$formula, 
            c(fit$y[which(censoring==1)] + rgl(nrow(fit$Fitted[which(censoring==1),,drop=F]), 
            fit$"Estimated parameters"[match(c("L1", "L2", "L3", 
                "L4"), names(fit$"Estimated parameters"))], param = fit$param),fit$y[which(censoring==0)] + rgl(nrow(fit$Fitted[which(censoring==0),,drop=F]), 
            fit$censor.gld.values, param = fit$param))~.), fit, censoring, data = na.omit(subset(data, 
            select = all.vars(fit$formula))), param = fit$param, 
            fun = fun, init = init,fit=fit,censor.type=censor.type,adj.int = adj.int, GLD.adj = GLD.adj, adj.censor = adj.censor)
    }
    fit.bc.f <- fun.simu.bias.correct.alt(fit.simu.f, fit)
    result <- list(fit, simu.result=fit.simu.f, simu.bias.correct.result=fit.bc.f)
    if (summary.plot == TRUE) {
        summaryGraphics.gld.surv.lm(result,range=range)
    }
    return(result)
}

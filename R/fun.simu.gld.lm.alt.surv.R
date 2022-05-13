fun.simu.gld.lm.alt.surv <-
function (n.simu, formula, fit.obj, censoring, data, param, maxit = 20000, fun,
    method = "Nelder-Mead", init = NULL,fit=NULL,censor.type="right",adj.int = FALSE, GLD.adj = FALSE, adj.censor = TRUE) 
{
  
    index.var <- 1:(match("L1", names(fit.obj[[3]])) - 1)
    pb <- txtProgressBar(min = 0, max = n.simu, style = 3)
    r <- lapply(1:n.simu, function(i, formula, fit.obj, index.var, 
        data, param, fun, maxit, init,fit,censor.type,adj.int,GLD.adj,adj.censor,censoring) {
        setTxtProgressBar(pb, i)
        GLD.lm.simu.surv(formula, censoring=censoring, init.coeff = c(fit.obj[[3]][index.var]), 
            init.resid = fit.obj[[3]][(max(index.var) + 1):(max(index.var) + 
                4)], data = data, param = param, fun = fun, method = method, 
            maxit = maxit, fit = fit.obj, 
            init = init,censor.type=censor.type,adj.int = adj.int, GLD.adj = GLD.adj, adj.censor = adj.censor)
    }, formula, fit.obj, index.var, data, param, fun, maxit, 
        init,fit,censor.type,adj.int,GLD.adj,adj.censor,censoring)
    close(pb)
    r <- do.call("rbind", r)
    return(r)
}

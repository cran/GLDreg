\name{GLD.lm.full.surv}
\alias{GLD.lm.full.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function fits a GLD Accelerated Failure Time regression linear model and 
conducts simulations to display the statistical properties of estimated 
coefficients
}
\description{
The function is an extension of \code{\link{GLD.lm.surv}} and defaults to
1000 simulation runs, coefficients and statistical properties of coefficients
can be plotted as part of the output.
}
\usage{
GLD.lm.full.surv(formula, censoring, data, param, maxit = 20000, fun, 
method = "Nelder-Mead", range = c(0.01, 0.99), n.simu = 1000, 
summary.plot = FALSE, init = NULL, alpha = 0.05, censor.type = "right", 
adj.int = FALSE, GLD.adj = FALSE, adj.censor = TRUE, keep.uncen = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
A symbolic expression of the model to be fitted, similar to the formula
argument in \code{\link{lm}}, see \code{\link{formula}} for more information
}
  \item{censoring}{
1=Event, 0= Censored
}
  \item{data}{
Dataset containing variables of the model
}
  \item{param}{
Can be "rs", "fmkl" or "fkml"
}
  \item{maxit}{
Maximum number of iterations for numerical optimisation 
}
  \item{fun}{
If param="fmkl" or "fkml", this can be one of \code{\link[GLDEX]{fun.RMFMKL.ml.m}}, 
\code{\link[GLDEX]{fun.RMFMKL.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml).

If param="rs", this can be one of \code{\link[GLDEX]{fun.RPRS.ml.m}}, 
\code{\link[GLDEX]{fun.RPRS.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml).
}
  \item{method}{
Defaults to "Nelder-Mead" algorithm, can also be "SANN" but this is a lot slower
and may not as good
}
 \item{range}{
The is the quantile range to plot the QQ plot, defaults to 0.01 and 0.99 to
avoid potential problems with extreme values of GLD which might be -Inf or Inf.
}
  \item{n.simu}{
Number of simulations, defaults to 1000.
}
  \item{summary.plot}{
If TRUE present graphical display of model fitted.
}
  \item{init}{
Initial values to start optimization process.
}
  \item{alpha}{
Significant level of goodness of fit test.
}
\item{censor.type}{
Can be " right" of "left censored.
}
  \item{adj.int}{
Adjust intercept in final output?
}
  \item{GLD.adj}{
Adjust GLD fitted to have theoretical zero mean?
}
  \item{adj.censor}{
Adjust censoring?
}
  \item{keep.uncen}{
Keep uncensored values?
}
}

\value{

 \item{Message}{Short description of estimation method used and whether the 
 result converged}
 \item{Bias Correction}{Bias correction used to ensure the line has zero mean 
 residuals}
 \item{Estimated parameters}{A set of estimate coefficients from GLD regression}
 \item{Fitted}{Predicted response value from model}
 \item{Residual}{Residual of model}
 \item{formula}{Formula used in the model}
 \item{param}{Specify whether RS/FKML/FMKL GLD was used}
 \item{y}{The response variable} 
 \item{x}{The explanatory variable(s)}
 \item{fun}{GLD fitting function used in the computation process, outputted for 
 internal programming use} 
 \item{censoring}{Censoring data} 
 \item{AIC.full}{AIC results} 
 \item{BIC.full}{BIC results} 
 \item{censor.gld.values}{Result of GLD fit, including censoring} 
 \item{simu.result}{Result of simulation for all coefficeints in the model} 
 \item{censor.gld.values}{Result of GLD fit, including censoring} 
 \item{simu.bias.correct.result}{Bias corrected simulation results}
}
\references{
Su (2021) "Flexible Parametric Accelerated Failure Time Model" Journal of 
Biopharmaceutical Statistics Volume 31, 2021 - Issue 5
}
\author{
Steve Su
}

\seealso{
\code{\link{GLD.lm.full}}, \code{\link{GLD.quantreg}}, \code{\link{GLD.lm}}, 
\code{\link{GLD.lm.surv}}
}

\examples{
\dontrun{

actg.rs<-GLD.lm.full.surv(log(time)~factor(txgrp)+hemophil+cd4+priorzdv+age,
censoring=actg[which(actg$txgrp!=3 & actg$txgrp!=4),]$censor, 
data=actg[which(actg$txgrp!=3 & actg$txgrp!=4),],
param="rs",fun=fun.RPRS.ml.m,summary.plot=F,n.simu=1000)

}}

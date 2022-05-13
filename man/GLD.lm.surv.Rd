\name{GLD.lm.surv}
\alias{GLD.lm.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function fits a GLD Accelerated Failure Time Model for Survival Data
}
\description{
Similar to \code{\link{GLD.lm}}, this function fits an Accelerated Failure Time 
Model using RS/FKML GLDs.
}
\usage{
GLD.lm.surv(formula, censoring, data, param, maxit = 20000, fun, 
            method = "Nelder-Mead", diagnostics = TRUE, range = c(0.01, 0.99), 
            init = NULL, alpha = 0.05, censor.type = "right", adj.int = TRUE, 
            GLD.adj = FALSE, adj.censor = TRUE, keep.uncen = TRUE)
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
%%     ~~Describe \code{maxit} here~~
}
  \item{fun}{
If param="fmkl" or "fkml", this can be one of \code{\link{fun.RMFMKL.ml.m}}, 
\code{\link{fun.RMFMKL.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml).

If param="rs", this can be one of \code{\link{fun.RPRS.ml.m}}, 
\code{\link{fun.RPRS.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml).
}
  \item{method}{
Defaults to "Nelder-Mead" algorithm, can also be "SANN" but this is a lot slower
and may not as good
}
  \item{diagnostics}{
Defaults to TRUE, which computes Kolmogorov-Smirnoff test, Kolmogorov-Smirnoff 
Resample test, Data drive smooth test and do QQ plot on non censored data.

}
  \item{range}{
The is the quantile range to plot the QQ plot, defaults to 0.01 and 0.99 to
avoid potential problems with extreme values of GLD which might be -Inf or Inf.

}
  \item{init}{
Choose a different set of initial values to start the optimisation process. This
can either be full set of parameters including GLD parameter estimates, or it
can just be the coefficient estimates of the regression model. 
}
  \item{alpha}{
Significant level of goodness of fit test.
}
  \item{censor.type}{
Can be" right" of "left censored.
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
\code{\link{GLD.lm.full.surv}}
}

\examples{

\dontrun{
library(mlr3proba)

actg320.rs<-GLD.lm.full(log(time)~factor(txgrp)+hemophil+cd4+priorzdv+age,
censoring=actg320[which(actg320$txgrp!=3 & actg320$txgrp!=4),]$censor, 
data=actg320[which(actg320$txgrp!=3 & actg320$txgrp!=4),],
param="rs",fun=fun.RPRS.ml.m)

}
}
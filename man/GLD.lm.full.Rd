\name{GLD.lm.full}
\alias{GLD.lm.full}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function fits a GLD regression linear model and conducts simulations to 
display the statistical properties of estimated coefficients
}
\description{
The function is an extension of \code{\link{GLD.lm}} and defaults to
1000 simulation runs, coefficients and statistical properties of coefficients
can be plotted as part of the output.
}
\usage{
GLD.lm.full(formula, data, param, maxit = 20000, fun, method = "Nelder-Mead", 
range = c(0.01, 0.99), n.simu = 1000, summary.plot = TRUE, init = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A symbolic expression of the model to be fitted, similar to the formula
argument in \code{\link{lm}}, see \code{\link{formula}} for more information
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
If param="fmkl" or "fkml", this can be one of \code{\link{fun.RMFMKL.ml.m}}, 
\code{\link{fun.RMFMKL.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml) and 
\code{\link{fun.RMFMKL.lm}} for L moment matching.

If param="rs", this can be one of \code{\link{fun.RPRS.ml.m}}, 
\code{\link{fun.RPRS.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml) and 
\code{\link{fun.RPRS.lm}} for L moment matching.

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
Number of times to repeat the simulation runs, defaults to 1000.
}
  \item{summary.plot}{
Whether to plot the coefficients graphically, defaults to TRUE.
}
 \item{init}{
Choose a different set of initial values to start the optimisation process. This
can either be full set of parameters including GLD parameter estimates, or it
can just be the coefficient estimates of the regression model. 
}
}
\details{
This function usually takes some time to run, as it involves refitting the 
GLD regression model many times, the progress of the simulation is outputted
to the R screen, so users can guage the progress of the computation.
}
\value{

\item{[[1]]}{Output of \code{\link{GLD.lm}}}
\item{[[2]]}{A matrix showing the bias adjustment, coefficents of the model,
parameters of GLD and whether the result converged at each run}
\item{[[3]]}{Adjusted simulation result so that the empirical mean of 
coefficients is the same as the estimated parameters obtained in 
\code{\link{GLD.lm}}}

}
\references{
Su (2015) "Flexible Parametric Quantile Regression Model" Statistics &
Computing  May 2015, Volume 25, Issue 3, pp 635-650
}
\author{
Steve Su
}
\seealso{
\code{\link{GLD.lm}}, \code{\link{GLD.quantreg}},  
\code{\link{summaryGraphics.gld.lm}}
}
\examples{

## Dummy example

## Create dataset

set.seed(10)

x<-rnorm(200,3,2)
y<-3*x+rnorm(200)

dat<-data.frame(y,x)

## Fit FKML GLD regression with 3 simulations

fit<-GLD.lm.full(y~x,data=dat,fun=fun.RMFMKL.ml.m,param="fkml",n.simu=3)

\dontrun{
## Extract the Engel dataset 

library(quantreg)
data(engel)

## Fit a full GLD regression

engel.fit.full<-GLD.lm.full(foodexp~income,data=engel,param="fmkl",
fun=fun.RMFMKL.ml.m)

## Extract the mammals dataset 
library(MASS)

## Fit a full GLD regression

mammals.fit.full<-GLD.lm.full(log(brain)~log(body),data=mammals,param="fmkl",
fun=fun.RMFMKL.ml.m)

## Using quantile regression coefficients as starting values
library(quantreg)

mammals.fit1.full<-GLD.lm.full(log(brain)~log(body),data=mammals,param="fmkl",
fun=fun.RMFMKL.ml.m, init=rq(log(brain)~log(body),data=mammals)$coeff)

## Using the result of mammals.fit.full as initial values

mammals.fit2.full<-GLD.lm.full(log(brain)~log(body),data=mammals,param="fmkl",
fun=fun.RMFMKL.ml.m, init=mammals.fit1.full[[1]][[3]])
}
}
\keyword{model}


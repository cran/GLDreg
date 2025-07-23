\name{GLD.lm}
\alias{GLD.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function fits a GLD regression linear model
}
\description{
Similar to \code{\link{lm}}, this function fits a linear model using RS/FKML
GLDs and assess the goodness of fit of GLD with respect to the data via qq plot
and Kolmogorov-Smirnoff (KS) test. Note that the use of KS test when parameters 
of a distribution are estimated from data is generally frowned upon. This is because 
one often gets inflated p-value with increased type II error due to the fact that the 
KS test requires independence between test sample and parameters of distribution.
Therefore, the the resample KS test over 1000 simulation runs from GLDEX package is 
probably a more reasonable measure. It is probably reasonable to consider the resample 
KS test may in fact decrease the p-values, as testing is done on resampled data from 
fitted distribution so there is a certain degree of inaccuracy there. The provision of 
these results is to give some indication of optimistic and pessimistic goodness
of fit measure, as currently, there is an absence of a specialised GLD goodness of fit test. 
A generic Data Driven Smooth Test from ddst library in R is also incorporated to assess goodness of fit.

When in doubt, QQ plot should always be considered ahead of these results.
}
\usage{
GLD.lm(formula, data, param, maxit = 20000, fun, method = "Nelder-Mead", 
diagnostics = TRUE, range = c(0.01, 0.99), init = NULL, alpha = 0.05)
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
If param="fmkl" or "fkml", this can be one of \code{\link[GLDEX]{fun.RMFMKL.ml.m}}, 
\code{\link[GLDEX]{fun.RMFMKL.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml) and 
\code{\link[GLDEX]{fun.RMFMKL.lm}} for L moment matching.

If param="rs", this can be one of \code{\link[GLDEX]{fun.RPRS.ml.m}}, 
\code{\link[GLDEX]{fun.RPRS.ml}}, for maximum
likelihood estimation (*.ml.m is a faster implementation of *.ml) and 
\code{\link[GLDEX]{fun.RPRS.lm}} for L moment matching.

}
  \item{method}{
Defaults to "Nelder-Mead" algorithm, can also be "SANN" but this is a lot slower
and may not as good
}
  \item{diagnostics}{
Defaults to TRUE, which computes Kolmogorov-Smirnoff test and do QQ plot 
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
Significant level of KS test.
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
}
\references{
Su (2015) "Flexible Parametric Quantile Regression Model" Statistics &
Computing  May 2015, Volume 25, Issue 3, pp 635-650
}
\author{
Steve Su
}

\seealso{
\code{\link{GLD.lm.full}}, \code{\link{GLD.quantreg}}
}
\examples{

## Dummy example

library(GLDEX)

## Create dataset

set.seed(10)

x<-rnorm(200,3,2)
y<-3*x+rnorm(200)

dat<-data.frame(y,x)

## Fit a FKML GLD regression

example<-GLD.lm(y~x,data=dat,fun=fun.RMFMKL.ml.m,param="fkml")

\dontrun{

## Extract the Engel dataset 
library(quantreg)
data(engel)

## Fit GLD Regression
engel.fit<-GLD.lm(foodexp~income,data=engel,param="fmkl",fun=fun.RMFMKL.ml.m)

## Extract the mammals dataset 
library(MASS)

mammals.fit<-GLD.lm(log(brain)~log(body),data=mammals,param="rs",
fun=fun.RPRS.lm)

## Using quantile regression coefficients as starting values
library(quantreg)

mammals.fit1<-GLD.lm(log(brain)~log(body),data=mammals,param="rs",
fun=fun.RPRS.lm,init=rq(log(brain)~log(body),data=mammals)$coeff)

# As an exercise, use the result from mammals.fit1 as initial values

GLD.lm(log(brain)~log(body),data=mammals,param="rs",
fun=fun.RPRS.lm,init=mammals.fit1[[3]])
}
}
\keyword{model}


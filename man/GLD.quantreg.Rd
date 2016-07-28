\name{GLD.quantreg}
\alias{GLD.quantreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a GLD quantile regression parametrically or non parametrically
}
\description{
The GLD quantile regression can be: 1) Fixed intercept, allowing all other
coefficients to vary, 2) Only intercept is allowed to vary and 3) All 
coefficients can vary. Minimisation is achieved numerically through least
squares between the proportion of estimated GLD error distribution below zero 
versus the specified quantile for parametric approach. For non parametric 
approach, minimisation is achieved using a least squares approach to find a
q-th quantile GLD line such that the percentage of observations below
the line corresponds to the q-th quantile. }
\usage{
GLD.quantreg(q, fit.obj, intercept = "", slope = "", emp=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
Specify the quantile (range 0 to 1) line
}
  \item{fit.obj}{
An object from \code{\link{GLD.lm.full}}
}
  \item{intercept}{
Can either be "fixed" or left blank, blank indicates this parameter is
allowed to vary in quantile line estimation
}
  \item{slope}{
Can either be \code{"fixed"} or left blank, blank indicates this parameter is
allowed to vary in quantile line estimation
}

  \item{emp}{
Can either be \code{TRUE} (non parametric GLD quantile regression) or 
\code{FALSE} (parametric GLD quantile regression), defaults to \code{FALSE}
}
}
\details{
This is a wrapper function for \code{\link{fun.gld.all.vary}},
\code{\link{fun.gld.slope.fixed.int.vary}}, 
\code{\link{fun.gld.slope.vary.int.fixed}}.
}
\value{
A matrix showing the estimated coefficients for the specified quantile
regression model, the objective function value and whether convergence is 
reached in the optimisation process. A value of 0 indicates convergence is 
reached. The convergence value is the same as the one from the \code{\link{optim}} 
function.
}
\references{
Su (2015) "Flexible Parametric Quantile Regression Model" Statistics &
Computing  May 2015, Volume 25, Issue 3, pp 635-650
}
\author{
Steve Su
}
\seealso{
\code{\link{GLD.lm.full}},\code{\link{fun.plot.q}},  
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

## Find median regression, use empirical method

med.fit<-GLD.quantreg(0.5,fit,slope="fixed",emp=TRUE)

\dontrun{

## Extract the Engel dataset 

library(quantreg)
data(engel)

## Fit GLD Regression along with simulations

engel.fit.all<-GLD.lm.full(foodexp~income,data=engel,
param="fmkl",fun=fun.RMFMKL.ml.m)

## Fit parametric GLD quantile regression from 0.1 to 0.9, with equal spacings 
## between quantiles

result<-GLD.quantreg(seq(0.1,.9,length=9),engel.fit.all,intercept="fixed")

## Non parametric quantile regression

GLD.quantreg(seq(0.1,.9,length=9),engel.fit.all,intercept="fixed",emp=T)

}
}

\keyword{model}









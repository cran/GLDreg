\name{fun.plot.q}
\alias{fun.plot.q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
2-D Plot for Quantile Regression lines
}
\description{
This function plots quantile regression lines from \code{\link{GLD.lm}} and
one of \code{\link{fun.gld.slope.vary.int.fixed}},
\code{\link{fun.gld.slope.fixed.int.vary}}, 
\code{\link{fun.gld.slope.fixed.int.vary.emp}}, 
\code{\link{fun.gld.all.vary.emp}}, \code{\link{fun.gld.all.vary}}, 
\code{\link{fun.gld.slope.vary.int.fixed.emp}}, \code{\link{GLD.quantreg}}.
}
\usage{
fun.plot.q(x, y, fit, quant.info, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numerical vector of explanatory variable
}
  \item{y}{
A numerical vector of response variable 
}
  \item{fit}{
An object from \code{\link{GLD.lm}}
}
  \item{quant.info}{
An object from one of \code{\link{fun.gld.slope.vary.int.fixed}},
\code{\link{fun.gld.slope.fixed.int.vary}}, 
\code{\link{fun.gld.slope.fixed.int.vary.emp}}, 
\code{\link{fun.gld.all.vary.emp}}, \code{\link{fun.gld.all.vary}}, 
\code{\link{fun.gld.slope.vary.int.fixed.emp}}, \code{\link{GLD.quantreg}}
}
  \item{\dots}{
Additional arguments to be passed to plot function, such as axis labels and 
title of the graph
}
}
\details{
This is intended to plot only two variables, for quantile regression involving
more than one explanatory variable, consider plotting the actual values versus 
fitted values by fitting a secondary GLD quantile model between actual and 
fitted values.
}
\value{
A graph showing quantile regression lines
}
\references{
Su (In Press) "Flexible Parametric Quantile Regression Model" Statistics &
Computing  
}
\author{
Steve Su
}

\examples{

## Dummy example

library(GLDEX)

## Create dataset

set.seed(10)

x<-rnorm(200,3,2)
y<-3*x+rnorm(200)

dat<-data.frame(y,x)

## Fit FKML GLD regression with 5 simulations

fit<-GLD.lm.full(y~x,data=dat,fun=fun.RMFMKL.ml.m,param="fkml",n.simu=5)

## Find median regression, use empirical method

med.fit<-GLD.quantreg(0.5,fit,slope="fixed",emp=TRUE)

fun.plot.q(x=x,y=y,fit=fit[[1]],med.fit, xlab="x",ylab="y")

\donttest{

## Plot result of quantile regression

## Extract the Engel dataset 

library(quantreg)
data(engel)

## Fit GLD Regression along with simulations

engel.fit.all<-GLD.lm.full(foodexp~income,data=engel,
param="fmkl",fun=fun.RMFMKL.ml.m)

## Fit quantile regression from 0.1 to 0.9, with equal spacings between 
## quantiles

result<-GLD.quantreg(seq(0.1,.9,length=9),engel.fit.all,intercept="fixed")

## Plot the quantile regression lines

fun.plot.q(x=engel$income,y=engel$foodexp,fit=engel.fit.all[[1]],result,
xlab="income",ylab="Food Expense")
}
}
\keyword{hplot}


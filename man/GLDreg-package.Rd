\name{GLDreg-package}
\alias{GLDreg-package}
\alias{GLDreg}
\docType{package}
\title{
This package fits standard and quantile regression models using RS and FMKL/FKML 
generalised lambda distributions via maximum likelihood estimation and L moment 
matching.
}
\description{
Owing to the rich shapes of GLDs, GLD standard/quantile regression is a 
competitive flexible model compared to standard/quantile regression.  The 
proposed method has some major advantages: 1) it provides a reference line which 
is very robust to outliers with the attractive property of zero mean residuals 
and 2) it gives a unified, elegant quantile regression model from the reference 
line with smooth regression coefficients across different quantiles. The 
goodness of fit of the proposed model can be assessed via QQ plots and the 
Kolmogorov-Smirnov test, to ensure the appropriateness of the statistical 
inference under consideration. Statistical distributions of coefficients of the
GLD regression line are obtained using simulation, and interval estimates are 
obtained directly from simulated data. 
}
\details{
\tabular{ll}{
Package: \tab GLDreg\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2014-03-04\cr
License: \tab GPL (>= 3.0)\cr
}
The primary fitting function for GLD regression model is
\code{\link{GLD.lm.full}}. The output of \code{\link{GLD.lm.full}} can then be
passed to \code{\link{summaryGraphics.gld.lm}} to display coefficients of GLD
regression model graphically. Once a GLD reference model is obtained,  
quantile regression is obtained using \code{\link{GLD.quantreg}}.
}
\author{
Steve Su <allegro.su@gmail.com>
}
\references{
Su (2014) "Flexible Parametric Quantile Regression Model" Statistics &
Computing  
}
\seealso{
\pkg{GLDEX}
}
\examples{

## Dummy example

## Create dataset

set.seed(10)

x<-rnorm(200,3,2)
y<-3*x+rnorm(200)

dat<-data.frame(y,x)

## Fit a FKML GLD regression

example<-GLD.lm(y~x,data=dat,fun=fun.RMFMKL.ml.m,param="fkml")

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

## Plot coefficient summary

summaryGraphics.gld.lm(engel.fit.all)

## Fit quantile regression from 0.1 to 0.9, with equal spacings between 
## quantiles

result<-GLD.quantreg(seq(0.1,.9,length=9),engel.fit.all,intercept="fixed")

## Plot quantile regression lines

fun.plot.q(x=engel$income,y=engel$foodexp,fit=engel.fit.all[[1]],result,
xlab="income",ylab="Food Expense")

## Sometimes the maximum likelihood estimation may fail, for example when 
## minimum/maximum support of GLD is exactly at the minimum/maximum value of the 
## dataset, if this the case, try to use the L-moment matching method.

engel.fit.all<-GLD.lm.full(foodexp~income,data=engel,
param="fmkl",fun=fun.RMFMKL.lm)

}
}

\keyword{model}
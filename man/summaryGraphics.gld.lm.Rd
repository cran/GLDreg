\name{summaryGraphics.gld.lm}
\alias{summaryGraphics.gld.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graphical display of output from \code{\link{GLD.lm.full}}
}
\description{
This function display the coefficients and the distribution of coefficients 
obtained from GLD regression model.
}
\usage{
summaryGraphics.gld.lm(overall.fit.obj, alpha = 0.05, label = NULL, 
ColourVersion = TRUE, diagnostics = TRUE, range = c(0.01, 0.99))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{overall.fit.obj}{
An object from \code{\link{GLD.lm.full}}
}
  \item{alpha}{
Specifying the range of interval for the coefficients, default is 0.05, which
specifies a 95\% interval. 
}
  \item{label}{
A character vector indicating the labelling for the coeffiients
}
  \item{ColourVersion}{
Whether to display colour or not, default is TRUE, if set as FALSE, a black
and white plot is given. This is only applicable to the coefficient summary
graph and has no effect on QQ plots.
}
  \item{diagnostics}{
If TRUE, then QQ plot will be given along with Kolmogorov-Smirnoff test results
}
  \item{range}{
The is the quantile range to plot the QQ plot, defaults to 0.01 and 0.99 to
avoid potential problems with extreme values of GLD which might be -Inf or Inf.
}
}
\details{
The reason QQ plots are not displayed in black and white even if 
ColourVersion is set to FALSE is because the colour is necessary in those plots
for clarity of display. 
}
\value{
Graphics displaying coefficient summary and diagnostic plot (if chosen)
}
\references{
Su (In Press) "Flexible Parametric Quantile Regression Model" Statistics &
Computing  
}
\author{
Steve Su
}
\seealso{
\code{\link{GLD.lm.full}}
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

summaryGraphics.gld.lm(fit,ColourVersion=FALSE,diagnostic=FALSE)

\donttest{
## Extract the Engel dataset 

library(quantreg)
data(engel)

## Fit a full GLD regression

engel.fit.full<-GLD.lm.full(foodexp~income,data=engel,param="fmkl",
fun=fun.RMFMKL.ml.m)

## Plot coefficient summary

summaryGraphics.gld.lm(engel.fit.full,ColourVersion=FALSE,diagnostic=FALSE)

summaryGraphics.gld.lm(engel.fit.full)

## Extract the mammals dataset 
library(MASS)

## Fit a full GLD regression

mammals.fit.full<-GLD.lm.full(log(brain)~log(body),data=mammals,param="fmkl",
fun=fun.RMFMKL.ml.m)

## Plot coefficient summary

summaryGraphics.gld.lm(mammals.fit.full,label=c("intercept","log of body weight"))

}
}
\keyword{hplot}


\name{qqgld.default}
\alias{qqgld.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QQ plot for GLD
}
\description{
This is an updated QQ plot function for GLD comparing fitted distribution
with empirical data
}
\usage{
qqgld.default(y, vals, param, ylim, main = "GLD Q-Q Plot", 
xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", 
plot.it = TRUE, datax = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A vector of empirical data observations
}
  \item{vals}{
A vector representing four parameters of GLD
}
  \item{param}{
Can be "rs", "fmkl" or "fkml"
}
  \item{ylim}{
A vector of two numerical values, specifying the upper and lower bound
of y axis
}
  \item{main}{
Title of the qq plot
}
  \item{xlab}{
Label for X axis
}
  \item{ylab}{
Label for Y axis
}
  \item{plot.it}{
Whether to plot the QQ plot, default is TRUE
}
  \item{datax}{
Whether data values should be on x axis, default is FALSE
}
  \item{\dots}{
Additional graphical parameters  
}
}
\details{
This is an adaptation of the default qq plot in R
}
\value{
A list with components:

\item{x}{	
The x coordinates of the points that were/would be plotted  }

\item{y}{	
The original y vector, i.e., the corresponding y coordinates including NAs.
}
 }
\author{
R, with modifications from Steve Su
}
\seealso{
\code{\link{qqplot.gld}}, \code{\link{qqplot.gld.bi}} 
}
\examples{
x<-rnorm(100)
fit1<-fun.RMFMKL.ml.m(x)
qqgld.default(x,fit1,param="fmkl")
}
\keyword{hplot}


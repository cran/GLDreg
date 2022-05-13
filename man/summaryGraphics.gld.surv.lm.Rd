\name{summaryGraphics.gld.surv.lm}
\alias{summaryGraphics.gld.surv.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graphical display of output from \code{\link{GLD.lm.full.surv}}
}
\description{
This function display the coefficients and the distribution of coefficients 
obtained from GLD Accelerated Failure Time regression model.}
\usage{
summaryGraphics.gld.surv.lm(overall.fit.obj, alpha = 0.05, label = NULL, 
                            ColourVersion = TRUE, diagnostics = TRUE, 
                            range = c(0.01, 0.99), exp = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{overall.fit.obj}{
An object from \code{\link{GLD.lm.full.surv}}
}
  \item{alpha}{
Specifying the range of interval for the coefficients, default is 0.05, which
specifies a 95\% interval. This also specifies the significance level of KS 
test.
}
  \item{label}{
A character vector indicating the labelling for the coefficients
}
  \item{ColourVersion}{
Whether to display colour or not, default is TRUE, if set as FALSE, a black
and white plot is given. This is only applicable to the coefficient summary
graph and has no effect on QQ plots.

}
  \item{diagnostics}{
If TRUE, then QQ plot will be given along with various goodness of fit test 
results}

  \item{range}{
The is the quantile range to plot the QQ plot, defaults to 0.01 and 0.99 to
avoid potential problems with extreme values of GLD which might be -Inf or Inf.
}
  \item{exp}{
If TRUE, Exponentiate the coefficients
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
Su (2021) "Flexible Parametric Accelerated Failure Time Model" Journal of 
Biopharmaceutical Statistics Volume 31, 2021 - Issue 5
}
\author{
Steve Su
}

\seealso{
\code{\link{GLD.lm.full.surv}}}
\examples{

\dontrun{

library(mlr3proba)

actg320.rs<-GLD.lm.full.surv(log(time)~factor(txgrp)+hemophil+cd4+priorzdv+age,
censoring=actg320[which(actg320$txgrp!=3 & actg320$txgrp!=4),]$censor, 
data=actg320[which(actg320$txgrp!=3 & actg320$txgrp!=4),],
param="rs",fun=fun.RPRS.ml.m,summary.plot=F,n.simu=1000)

summaryGraphics.gld.surv.lm(actg320.rs,label=c("(Intercept)",
"IDV versus no IDV","Hemophiliac","Baseline CD4",
"Months of prior \n ZDV use","Age"),exp="TRUE")

}
}

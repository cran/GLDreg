\name{fun.mean.convert}
\alias{fun.mean.convert}

\title{
Convert a RS or FKML GLD into RS or FKML GLD to the desired theoretical mean by
changing only the first parameter
}
\description{
A simple transformation of altering the location of RS/FKML GLD so that
the theoretical mean is altered to the level specified. Only the first parameter
of RS/FKML GLD is altered.
}
\usage{
fun.mean.convert(x, param, val = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of four values representing Lambda 1, Lambda 2, Lambda 3 and Lambda 4
of RS/FKML GLD.
}
  \item{param}{
Can be "rs" or "fmkl" or "fkml"
}
  \item{val}{
The targeted theoretical mean
}
}

\value{
A vector of four values representing Lambda 1, Lambda 2, Lambda 3 and Lambda 4
of the transformed RS/FKML GLD 
}

\author{
Steve Su
}
\note{
If finite first moment does not exist, original input values will be returned
}
\examples{

# Transform RS GLD with parameters 3,2,1,1 to mean of 0
fun.mean.convert(c(3,2,1,1),param="rs")

# Check that the desired outcome is achieved
fun.theo.mv.gld(0,2,1,1,param="rs")

# Transform RS GLD with parameters 3,2,1,1 to mean of 5
fun.mean.convert(c(3,2,1,1),param="fkml",5)

# Check that the desired outcome is achieved
fun.theo.mv.gld(5,2,1,1,param="fkml")

}
\keyword{univar}


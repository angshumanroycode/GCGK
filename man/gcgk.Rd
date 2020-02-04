\name{gcgk}
\alias{gcgk}
\title{
Dependency measure based on checkerboard copula and Gaussian kernel
}
\description{
Dependency measure that measures amount of association among several univariate
ordinal random variables.
}
\usage{
gcgk(X, parameters = NULL, B.parameters = 1000, 
normalize = c("unnormalized", "normalized"), preset = NULL)
}
\arguments{
  \item{X}{
A matrix whose rows represent observations.  
}
  \item{parameters}{
A numeric vector of parameters for the Gaussian kernel. If \code{NULL}, \code{gcgk} computes parameters automatically.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
  \item{normalize}{
It can either take \code{"normalized"} or \code{"unnormalized"} as input. Depending upon the input, it produces the results.
}
  \item{preset}{
If the list of numeric vectors is precalculated by the \code{gcgk.preset} function, supply that list here.
}
}
\details{
See reference.
}
\value{
A list containing numeric vectors of parameters, normalized or unnormalized gcgk values for corresponding parameters.
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
set.seed(0)
x <- rnorm(10)
y <- rep(rnorm(5),2)
z <- x+y+rnorm(10,sd=0.1)
X <- cbind(x,y,z)
gcgk(X)
}
\keyword{nonparametric}

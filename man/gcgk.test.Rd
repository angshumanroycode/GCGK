\name{gcgk.test}
\alias{gcgk.test}
\title{
Test of independence based on checkerboard copula and Gaussian kernel
}
\description{
Casrries out a test for independence among several univariate orinal random variables based on checkerboard copula and Gaussian kernel.
}
\usage{
gcgk.test(X, method = c("sum", "max", "fdr", "nsum", "nmax"),
                    B = 1000, parameters = NULL, B.parameters = 1000,
                    alpha = 0.05, preset = NULL)
}
\arguments{
  \item{X}{
A matrix whose rows represent observations.
}
  \item{method}{
A vector of strings containting any combination of these three strings: \code{"sum"}, \code{"max"}, \code{"fdr"}, \code{"nsum"} or \code{"nmax"}. Depending upon those combination, \code{gcgk.test} will perform those methods. For more details, see reference.
}
  \item{B}{
An integer value specifying the number of Monte-Carlo iterations to be used for determining the cutoff of the test.
}
  \item{parameters}{
If \code{NULL}, then \code{gcgk.test} automatically calculates the parameters. Otherwise it calculates preset using the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for determining the parameters.
}
  \item{alpha}{
A numeric value between 0 and 1 specifying the level of the test.
}
  \item{preset}{
If the list of numeric vectors is precalculated by the \code{gcgk.preset} function, supply that list here.
}
}
\details{
See reference.
}
\value{
A list depending upon the methods chosen.
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
gcgk.test(X)
}
\keyword{nonparametric}
\keyword{htest}

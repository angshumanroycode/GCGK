\name{gcgk.preset}
\alias{gcgk.preset}
\title{
A list of numeric vectors pre-calculated for \code{gcgk} or \code{gcgk.test} function
}
\description{
When \code{gcgk} or \code{gcgk.test} function has to be calculated repeatedly, we can save a lot of computation by pre-calculating some numeric vectors that \code{gcgk} or \code{gcgk.test} function uses by this \code{gcgk.preset} function.
}
\usage{
gcgk.preset(d, parameters = NULL, B.parameters = 1000, normalize = "normalized")
}
\arguments{
  \item{d}{
Number of univariate random variables.
}
  \item{parameters}{
If \code{NULL}, then \code{gcgk.preset} automatically calculates the parameters. Otherwise it calculates preset using the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
  \item{normalize}{
It can either take \code{"normalized"} or \code{"unnormalized"} as input. Depending upon the input, it produces list of vectors that is needed for \code{gcgk} or \code{gcgk.test} function.
}
}
\details{
See reference.
}
\value{
A list of numeric vectors.
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
#automated parameters
set.seed(0)
gcgk.preset(2,B.parameters=1000)

#manual parameters
gcgk.preset(2,parameters=c(0.3,0.5,1))
}
\keyword{nonparametric}

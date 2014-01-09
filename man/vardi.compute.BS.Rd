\name{vardi.compute.BS}
\alias{vardi.compute.BS}
\title{Compute B and S matrices in algorithm of Vardi (1996)}
\usage{
vardi.compute.BS(A, Y)
}
\arguments{
  \item{A}{routing matrix (m x k)}

  \item{Y}{matrix of link loads over time (m x n, one
  column per time)}
}
\value{
list containing two entries for the B and S matrices,
respectively
}
\description{
Function to compute B and S matrices for moment equations
of Vardi's method (1996). It's not particularly efficient,
but it works.
}
\references{
Y. Vardi. Network tomography: estimating source-destination
traffic intensities from link data. Journal of the American
Statistical Association, 91:365-377, 1996.
}
\seealso{
Other vardi: \code{\link{vardi.algorithm}};
\code{\link{vardi.iteration}}
}
\keyword{models}
\keyword{multivariate}


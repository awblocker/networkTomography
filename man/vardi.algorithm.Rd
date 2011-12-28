\name{vardi.algorithm}
\alias{vardi.algorithm}
\title{Run algorithm of Vardi (1996) given B and S matrices}
\usage{
  vardi.algorithm(A, Y, lambda, tol = 0.001)
}
\arguments{
  \item{A}{routing matrix (m x k)}

  \item{Y}{matrix of link loads over time (m x n, one
  column per time)}

  \item{lambda}{numeric vector of starting values for OD
  flows (length k)}

  \item{tol}{numeric tolerance for halting iterations}
}
\value{
  numeric vector of length k with estimated OD flows
}
\description{
  Runs moment-matching algorithm of Vardi (1996) until
  convergence
}
\references{
  Y. Vardi. Network tomography: estimating
  source-destination traffic intensities from link data.
  Journal of the American Statistical Association,
  91:365-377, 1996.
}
\seealso{
  Other vardi: \code{\link{vardi.compute.BS}},
  \code{\link{vardi.iteration}}
}
\keyword{models}
\keyword{multivariate}


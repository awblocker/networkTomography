\name{vardi.iteration}
\alias{vardi.iteration}
\title{Execute single iteration for algorithm of Vardi (1996)}
\usage{
  vardi.iteration(A, yBar, lambda, B, S)
}
\arguments{
  \item{A}{routing matrix (m x k)}

  \item{yBar}{numeric vector of mean link loads (length m)}
}
\value{
  list containing two entries for the B and S matrices,
  respectively
}
\description{
  Function to compute B and S matrices for moment equations
  of Vardi's method (1996). It's not particularly
  efficient, but it works.
}
\references{
  Y. Vardi. Network tomography: estimating
  source-destination traffic intensities from link data.
  Journal of the American Statistical Association,
  91:365-377, 1996.
}
\seealso{
  Other vardi: \code{\link{vardi.algorithm}},
  \code{\link{vardi.compute.BS}}
}
\keyword{models}
\keyword{multivariate}


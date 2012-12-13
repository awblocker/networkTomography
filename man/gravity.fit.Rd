\name{gravity.fit}
\alias{gravity.fit}
\title{Gravity estimation for a single time point}
\usage{
  gravity.fit(yt, srcDstInd)
}
\arguments{
  \item{yt}{length-m numeric vector of observed aggregate
  flows at time t}

  \item{srcDstInd}{list of source and destination flow
  indices corresponding to each point-to-point flow, as
  produced by \code{\link{getSrcDstIndices}}}
}
\value{
  xhat, a numeric vector of length k providing gravity
  estimates of the point-to-point flows of interest
}
\description{
  Gravity estimation for a single time point
}
\examples{
data(cmu)
srcDstInd <- getSrcDstIndices(cmu$A.full)
estimate <- gravity.fit(yt=cmu$Y.full[1,], srcDstInd=srcDstInd)
}
\seealso{
  Other gravity: \code{\link{gravity}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


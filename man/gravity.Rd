\name{gravity}
\alias{gravity}
\title{Run tomogravity estimation on complete time series of aggregate flows}
\usage{
gravity(Y, srcDstInd)
}
\arguments{
  \item{Y}{n x m matrix contain one vector of observed
  aggregate flows per row}

  \item{srcDstInd}{list of source and destination flow
  indices corresponding to each point-to-point flow, as
  produced by \code{\link{getSrcDstIndices}}}
}
\value{
Xhat, a n x k matrix containing a vector of estimated
point-to-point flows (for each time point) per row
}
\description{
Run tomogravity estimation on complete time series of
aggregate flows
}
\examples{
data(cmu)
srcDstInd <- getSrcDstIndices(cmu$A.full)
estimate <- gravity(Y=cmu$Y[1:3,], srcDstInd=srcDstInd)
}
\seealso{
Other gravity: \code{\link{gravity.fit}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


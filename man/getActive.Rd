\name{getActive}
\alias{getActive}
\title{Check for deterministically-known OD flows at single time}
\usage{
  getActive(y, A)
}
\arguments{
  \item{y}{numeric vector of link loads, dimension m}

  \item{A}{routing matrix of dimension m x k}
}
\value{
  logical vector of length k; TRUE for unknown OD flows,
  FALSE for known
}
\description{
  Uses xranges from limSolve to find
  deterministically-known OD flows
}
\examples{
data(onerouter)
getActive(Y[1,], A)
}
\keyword{algebra}


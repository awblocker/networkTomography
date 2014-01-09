\name{calcN}
\alias{calcN}
\title{Compute total traffic from a particular time.}
\usage{
calcN(yt, A1)
}
\arguments{
  \item{yt}{length-m numeric vectors of observed aggregate
  flows at a particular time}

  \item{A1}{m x m matrix containing the full-rank portion
  of the network's routing matrix, as supplied by
  \code{\link{decomposeA}}}
}
\description{
Compute total traffic from a particular time.
}
\examples{
data(bell.labs)
A.decomp <- decomposeA(bell.labs$A)
total.traffic <- calcN(yt=bell.labs$Y[1,], A1=A.decomp$A1)
total.traffic == sum(bell.labs$X[1,])
}
\keyword{algebra}


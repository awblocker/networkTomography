\name{getSrcDstIndices}
\alias{getSrcDstIndices}
\title{Find indices of source and destination for each point-to-point flow}
\usage{
  getSrcDstIndices(A)
}
\arguments{
  \item{A}{routing matrix of dimension m x k. This should
  be the reduced-rank version including all aggregate
  source and destination flows.}
}
\value{
  list consisting of two component, src and dst, which are
  integer vectors of length k containing the index (in y =
  A x) of the source and destination flows that each
  point-to-point flow is part of.
}
\description{
  Find indices of source and destination for each
  point-to-point flow
}
\examples{
data(cmu)
src.dst.ind <- getSrcDstIndices(cmu$A.full)
}
\keyword{algebra}


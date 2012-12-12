\name{decomposeA}
\alias{decomposeA}
\title{Compute pivoted decomposition of routing matrix A into full-rank and
remainder, as in Cao et al. 2000, via the QR decomposition.}
\usage{
  decomposeA(A)
}
\arguments{
  \item{A}{routing matrix of dimension m x k}
}
\value{
  list containing two matrices: A1 (m x m), a full-rank
  subset of the columns of A, and A2 (m x k - m), the
  remaining columns
}
\description{
  Compute pivoted decomposition of routing matrix A into
  full-rank and remainder, as in Cao et al. 2000, via the
  QR decomposition.
}
\keyword{algebra}


\name{diag_ind}
\alias{diag_ind}
\title{Make vector of 1-dimensional diagonal indices for square matrix}
\usage{
diag_ind(n)
}
\arguments{
  \item{n}{integer dimension of (square) matrix}
}
\value{
integer vector of length n with indices (unidimensional) of
square matrix
}
\description{
Compute vector of indices for efficient access to diagonal
of a square matrix
}
\examples{
ind <- diag_ind(5)
diag_mat(seq(5))[ind]
}
\seealso{
\code{\link{diag_mat}}
}
\keyword{array}


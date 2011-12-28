\name{diag_mat}
\alias{diag_mat}
\title{Make diagonal matrix from vector}
\usage{
  diag_mat(x)
}
\arguments{
  \item{x}{numeric vector for diagonal}
}
\value{
  matrix of size length(x) x length(x) with x along
  diagonal
}
\description{
  Build matrix with supplied vector on diagonal; this is
  much faster than diag due to the use of matrix instead of
  array
}
\examples{
diag_mat(seq(5))
}
\seealso{
  \code{\link{diag_ind}}
}
\keyword{array}


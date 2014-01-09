\name{agg}
\alias{agg}
\title{Function to aggregate results from matrix to matrix}
\usage{
agg(mat, q = c(0.05, 0.16, 0.5, 0.84, 0.95))
}
\arguments{
  \item{mat}{input numeric matrix to summarize}

  \item{q}{quantiles of mat's columns to provide in summary
  matrix}
}
\value{
matrix with each row corresponding to a summary measure and
each column corresponding to a column of mat
}
\description{
Defaults to mean, SD, limits, and given quantiles. Used to
limit memory consumption from MCMC runs.
}
\examples{
mat <- matrix(rnorm(5e3), ncol=5)
agg(mat)
}
\keyword{arith}
\keyword{manip}


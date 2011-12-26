\name{thin}
\alias{thin}
\title{Thinning vector of indices for MCMC}
\usage{
  thin(m, interval = 10)
}
\arguments{
  \item{m}{integer length of results}

  \item{interval}{thinning interval}
}
\value{
  integer vector of indices for thinning
}
\description{
  Returns a vector of indices with a given spacing for
  thinning MCMC results
}
\keyword{manip}
\keyword{ts}


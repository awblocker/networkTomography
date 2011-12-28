\docType{data}
\name{Y}
\alias{Y}
\title{1router link data from Cao et al. (2000)}
\description{
  Data from 4-node network with star topology collected
  from Bell Labs; used in Cao et al. (2000). The Y matrix
  contains only the link loads, formatted for analysis. The
  last (8th) link load has been removed to ensure that the
  resulting routing matrix is of full row rank.
}
\section{Structure}{
  The columns of the Y matrix correspond to individual link
  loads, and the rows correspond to observations. In this
  data, we have \code{A \%*\% t(X) == t(Y)}.
}
\references{
  J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
  Time-varying network tomography: router link data.
  Journal of the American Statistical Association,
  95:1063-75, 2000.
}
\seealso{
  \code{\link{onerouter}}, \code{\link{X}},
  \code{\link{tvec}}, \code{\link{A}}
}
\keyword{datasets}


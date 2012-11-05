\docType{data}
\name{cmu}
\alias{cmu}
\title{CMU data from Blocker & Airoldi (2011)}
\description{
  Data from the 12 node CMU network used in Blocker &
  Airoldi (2011). The OD flows are actual, observed traffic
  from a CMU network. The topology does not, however,
  correspond to the original network due to security
  considerations.
}
\section{Objects}{
  The list cmu, which contains several objects: \itemize{
  \item \code{A}, the routing matrix for this network
  (truncated for full row rank) \item \code{X}, a matrix of
  origin-destination flows formatted for analysis \item
  \code{Y}, a matrix of link loads formatted for analysis
  \item \code{A.full}, the routing matrix for this network
  without truncatation for full row rank) \item
  \code{Y.full}, a matrix of link loads corresponding to
  code{A.full} } In this data, we have \code{A \%*\% t(X)
  == t(Y)} and \code{A.full \%*\% t(X) == t(Y.full)}
}

\section{Variables}{
  The list cmu contains the following: \itemize{ \item The
  routing matrix \code{A}. The columns of this matrix
  correspond to individual OD flows (the columns of X), and
  its rows correspond to individual link loads (the columns
  of Y).  \item The OD matrix X. Columns correspond to
  individual OD flows, and the rows correspond to
  observations.  \item The link load matrix Y. Columns of
  the Y matrix correspond to individual link loads, and the
  rows correspond to observations.  \item The routing
  matrix \code{A.full}. This is the complete routing matrix
  before reduction for full row-rank.  \item The link load
  matrix Y.full, corresponding to A.full. }
}
\references{
  A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
  time series on a graph. Proceedings of the Twenty-Seventh
  Conference Annual Conference on Uncertainty in Artificial
  Intelligence (UAI-11) 51-60, 2011.
}
\keyword{datasets}


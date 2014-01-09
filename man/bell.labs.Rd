\docType{data}
\name{bell.labs}
\alias{bell.labs}
\title{Bell Labs 1router data from Cao et al. (2000)}
\usage{
bell.labs
}
\description{
Data from 4-node network with star topology collected from
Bell Labs; used in Cao et al. (2000).
}
\section{Objects}{
  The list bell.labs, which contains several objects:
  \itemize{ \item \code{A}, the routing matrix for this
  network (truncated for full row rank) \item \code{df}, a
  data.frame with all data \item \code{X}, a matrix of
  origin-destination flows formatted for analysis \item
  \code{Y}, a matrix of link loads formatted for analysis
  \item \code{tvec}, a vector of times } In this data, we
  have \code{A \%*\% t(X) == t(Y)}.
}

\section{Variables}{
  The list bell.labs contains the following: \itemize{
  \item The routing matrix \code{A}. The columns of this
  matrix correspond to individual OD flows (the columns of
  X), and its rows correspond to individual link loads (the
  columns of Y).  \item The data.frame \code{df},
  containing \itemize{ \item value, level of traffic
  recorded \item nme, name of flow or load \item method,
  whether flow was directly observered or inferred (all
  observed) \item time, time of observation \item od, flag
  for origin-destination vs. link loads \item orig, origin
  of flow or load \item dest, destination of flow or load
  \item node, node involved in flow or load } \item The OD
  matrix X. Columns correspond to individual OD flows, and
  the rows correspond to observations.  \item The link load
  matrix Y. Columns of the Y matrix correspond to
  individual link loads, and the rows correspond to
  observations.  \item The vector tvec, containing the time
  in decimal hours since midnight for each observation. }
}
\references{
J. Cao, D. Davis, S. Van Der Viel, and B. Yu. Time-varying
network tomography: router link data. Journal of the
American Statistical Association, 95:1063-75, 2000.
}
\keyword{datasets}


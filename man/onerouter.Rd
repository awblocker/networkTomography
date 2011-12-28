\docType{data}
\name{onerouter}
\alias{onerouter}
\title{1router data from Cao et al. (2000)}
\description{
  Data from 4-node network with star topology collected
  from Bell Labs; used in Cao et al. (2000).
}
\section{Objects}{
  This imports several objects: \itemize{ \item
  \code{onerouter}, a data.frame with all data \item
  \code{X}, a matrix of origin-destination flows formatted
  for analysis \item \code{Y}, a matrix of link loads
  formatted for analysis \item \code{tvec}, a vector of
  times \item \code{A}, the routing matrix for this network
  (truncated for full row rank) } In this data, we have
  \code{A \%*\% t(X) == t(Y)}.
}

\section{Variables}{
  The data.frame \code{onerouter} contains the following:
  \itemize{ \item value, level of traffic recorded \item
  nme, name of flow or load \item method, whether flow was
  directly observered or inferred (all observed) \item
  time, time of observation \item od, flag for
  origin-destination vs. link loads \item orig, origin of
  flow or load \item dest, destination of flow or load
  \item node, node involved in flow or load }
}
\references{
  J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
  Time-varying network tomography: router link data.
  Journal of the American Statistical Association,
  95:1063–75, 2000.
}
\seealso{
  \code{\link{X}}, \code{\link{Y}}, \code{\link{tvec}},
  \code{\link{A}}
}
\keyword{datasets}


\docType{data}
\name{abilene}
\alias{abilene}
\title{Abilene data from Fang et al. (2007)}
\usage{
abilene
}
\description{
Data from the 12 node Abilene network from Fang et al.
(2007). Both the OD flows and the topology correspond to
the actual network. This is the X1 dataset from the given
paper.
}
\section{Objects}{
  The list abilene, which contains several objects:
  \itemize{ \item \code{A}, the routing matrix for this
  network (truncated for full row rank) \item \code{X}, a
  matrix of origin-destination flows formatted for analysis
  \item \code{Y}, a matrix of link loads formatted for
  analysis \item \code{A.full}, the routing matrix for this
  network without truncatation for full row rank) \item
  \code{Y.full}, a matrix of link loads corresponding to
  code{A.full} } In this data, we have \code{A \%*\% t(X)
  == t(Y)} and \code{A.full \%*\% t(X) == t(Y.full)}
}

\section{Variables}{
  The list abilene contains the following: \itemize{ \item
  The routing matrix \code{A}. The columns of this matrix
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
J. Fang, Y. Vardi, and C.-H. Zhang. An iterative
tomogravity algorithm for the estimation of network
traffic. In R. Liu, W. Strawderman, and C.-H. Zhang,
editors, Complex Datasets and Inverse Problems: Tomography,
Networks and Beyond, volume 54 of Lecture Notes-Monograph
Series. IMS, 2007.
}
\keyword{datasets}


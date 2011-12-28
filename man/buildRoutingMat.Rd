\name{buildRoutingMat}
\alias{buildRoutingMat}
\title{Build routing matrices for linked star topologies; that is, a set of
star-topology networks with links between a subset of routers}
\usage{
  buildRoutingMat(nVec, Cmat)
}
\arguments{
  \item{nVec}{integer vector containing number of nodes in
  each sub-network (length m)}

  \item{Cmat}{matrix (m x m) containing a one for each
  linked sub-network; only upper triangular part is used}
}
\value{
  routing matrix of dimension at least 2*sum(nVec) x
  sum(nVec^2)
}
\description{
  Build routing matrices for linked star topologies; that
  is, a set of star-topology networks with links between a
  subset of routers
}
\examples{
nVec <- c(3, 3, 3)
Cmat <- diag(3)
Cmat[1,2] <- Cmat[2,3] <- 1
buildRoutingMat(nVec, Cmat)
}
\seealso{
  \code{\link{buildStarMat}}, which this function depends
  upon
}


\name{dobj.dxt.tomogravity}
\alias{dobj.dxt.tomogravity}
\title{Analytic gradient of objective function of Zhang et al. 2003}
\usage{
dobj.dxt.tomogravity(xt, yt, A, srcDstInd, lambda)
}
\arguments{
  \item{xt}{length-k numeric vector of point-to-point
  flows}

  \item{yt}{length-m numeric vector of observed aggregate
  flows}

  \item{A}{m x k routing matrix, yt = A xt}

  \item{srcDstInd}{list of source and destination flow
  indices corresponding to each point-to-point flow, as
  produced by \code{\link{getSrcDstIndices}}}

  \item{lambda}{regularization parameter for mutual
  information prior. Note that this is scaled by the
  squared total traffic in the objective function before
  scaling the mututal information prior.}
}
\value{
numeric vector of length k containing gradient of objective
function with respect to xt
}
\description{
Requires bounded optimization to maintain positive OD
flows, and only those flows that are not deterministically
zero should be included in the estimation.
}


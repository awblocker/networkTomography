\name{tomogravity.fit}
\alias{tomogravity.fit}
\title{Tomogravity estimation for a single time point using L-BFGS-B}
\usage{
tomogravity.fit(yt, A, srcDstInd, lambda, N = 1, normalize = FALSE,
  lower = 0, control = list())
}
\arguments{
  \item{yt}{length-m numeric vector of observed aggregate
  flows at time t}

  \item{A}{m x k routing matrix}

  \item{srcDstInd}{list of source and destination flow
  indices corresponding to each point-to-point flow, as
  produced by \code{\link{getSrcDstIndices}}}

  \item{lambda}{regularization parameter for mutual
  information prior. Note that this is scaled by the
  squared total traffic in the objective function before
  scaling the mututal information prior.}

  \item{N}{total traffic for normalization. Unused if
  normalized is FALSE.}

  \item{normalize}{If TRUE, xt and yt are scaled by N.
  Typically used in conjunction with calcN to normalize
  traffic to proportions, easing the tuning of lambda.}

  \item{lower}{Component-wise lower bound for xt in
  L-BFGS-B optimization.}

  \item{control}{List of control information for optim.}
}
\value{
A list as returned by optim, with element \code{par}
containing the estimated point-to-point flows and elementer
\code{gr} containing the analytic gradient evaluated at the
estimate.
}
\description{
Tomogravity estimation for a single time point using
L-BFGS-B
}
\examples{
data(cmu)
srcDstInd <- getSrcDstIndices(cmu$A.full)
estimate <- tomogravity.fit(yt=cmu$Y.full[1, ], A=cmu$A.full,
     srcDstInd=srcDstInd, lambda=0.01)
}
\seealso{
Other tomogravity: \code{\link{tomogravity}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


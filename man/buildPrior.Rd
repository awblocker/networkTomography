\name{buildPrior}
\alias{buildPrior}
\title{Construct prior from calibration model estimates}
\usage{
buildPrior(xHat, varHat, phiHat, Y, A, rho = 0.9, phiPriorDf = ncol(A)/2,
  backward = FALSE, lambdaMin = 1, ipfp.maxit = 1e+06, ipfp.tol = 1e-06)
}
\arguments{
  \item{xHat}{matrix (n x k) of estimates for OD flows from
  calibration model, one time point per row}

  \item{varHat}{matrix (n x k) of estimated variances for
  OD flows from calibration, one time point per row}

  \item{phiHat}{numeric vector (length n) of estimates for
  phi from calibration model}

  \item{Y}{matrix (n x l) of observed link loads, one time
  point per row}

  \item{A}{routing matrix (l x k) for network; must be of
  full row rank}

  \item{phiPriorDf}{numeric prior convolution parameter for
  independent inverse-gamma priors on phi_t}

  \item{rho}{numeric fixed autoregressive parameter for
  dynamics on lambda; see reference for details}

  \item{backward}{logical to activate construction of
  reversed prior (for smoothing applications)}

  \item{lambdaMin}{numeric value at which to floor
  estimated OD flows for prior construction}

  \item{ipfp.maxit}{integer maximum number of iterations
  for IPFP}

  \item{ipfp.tol}{numeric tolerance for convergence of IPFP
  iterations}
}
\value{
list containing priors for lambda and phi, consisting of:
\itemize{ \item mu, a matrix (n x k) containing the prior
means for the log-change in each lambda at each time \item
sigma, a matrix (n x k) containing the prior standard
deviations for the log-change in each lambda at each time
\item a list phi, containing the numeric prior \code{df}
and a numeric vector \code{scale} of length n }
}
\description{
Builds prior from appropriately structured output of the
calibration model from Blocker & Airoldi (2011). Handles
all formatting so result can be fed directly to
\code{\link{bayesianDynamicFilter}}.
}
\references{
A.W. Blocker and E.M. Airoldi. Deconvolution of mixing time
series on a graph. Proceedings of the Twenty-Seventh
Conference Annual Conference on Uncertainty in Artificial
Intelligence (UAI-11) 51-60, 2011.
}
\seealso{
Other bayesianDynamicModel:
\code{\link{bayesianDynamicFilter}};
\code{\link{move_step}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


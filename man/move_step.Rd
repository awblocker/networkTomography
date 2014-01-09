\name{move_step}
\alias{move_step}
\title{Move step of sample-resample-move algorithm for multilevel state-space model}
\usage{
move_step(y, X, tme, lambda, phi, lambdatm1, phitm1, prior, A, A1_inv, A2, rho,
  tau, m = ncol(X), l = nrow(A1_inv), k = length(lambda), ndraws = 10,
  minAccepts = 0, verbose = FALSE)
}
\arguments{
  \item{y}{numeric vector (length l) of observed link
  loads}

  \item{X}{matrix (m x k) of particles for OD flows, one
  particle per row, in pivoted order}

  \item{tme}{integer time index currently used in
  estimation}

  \item{lambda}{matrix (m x k) of particles for OD
  intensities, one particle per row, in pivoted order}

  \item{phi}{numeric vector (length m) of particles for
  phi}

  \item{lambdatm1}{lambda matrix (m x k) of particles for
  OD intensities from previous time, one particle per row,
  in pivoted order}

  \item{phitm1}{numeric vector (length m) of particles for
  phi from previous time}

  \item{prior}{list containing priors for hyperparameters;
  see \code{\link{bayesianDynamicFilter}} for details}

  \item{A}{routing matrix (l x k) for network}

  \item{A1_inv}{inverse of full-rank portion of routing
  matrix (l x l)}

  \item{A2}{remainder of routing matrix (l x k-l)}

  \item{rho}{numeric fixed autoregressive parameter for
  dynamics on lambda; see reference for details}

  \item{tau}{numeric fixed power parameter for variance
  structure on truncated normal noise; see reference for
  details}

  \item{m}{integer number of particles}

  \item{l}{integer number of observed link loads}

  \item{k}{integer number of OD flows to infer}

  \item{ndraws}{integer number of draws to perform (can be
  overriden by minAccepts)}

  \item{minAccepts}{integer minimum number of acceptances
  before results are returned; activates alternative
  stopping rule if >= 1}

  \item{verbose}{logical activates verbose diagnostic
  output}
}
\value{
list containing updated values of X, lambda, and phi
}
\description{
Function to execute single MCMC-based move step for
\code{\link{bayesianDynamicFilter}}. This can use two types
of stopping rules: number of iterations or number of
accepted moves for the X particles. The former is used by
default, but the latter adapts better to low acceptance
rates (sometimes with substantial computational cost). Most
updates in this algorithm are Metropolis-Hastings with
customized proposals.
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
\code{\link{buildPrior}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


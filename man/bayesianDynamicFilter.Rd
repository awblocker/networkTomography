\name{bayesianDynamicFilter}
\alias{bayesianDynamicFilter}
\title{Function for inference with multilevel state-space model}
\usage{
bayesianDynamicFilter(Y, A, prior, lambda0, sigma0, phi0, rho = 0.1,
  tau = 2, m = 1000, verbose = FALSE, Xdraws = 5 * m, Xburnin = m,
  Movedraws = 10, nThresh = 10, aggregate = FALSE, backward = FALSE,
  tStart = 1)
}
\arguments{
  \item{Y}{matrix (n x l) of observed link loads over time,
  one observation per row}

  \item{A}{routing matrix (l x k) for network; must be of
  full row rank}

  \item{prior}{list containing priors for lambda and phi;
  must have \itemize{ \item mu, a matrix (n x k) containing
  the prior means for the log-change in each lambda at each
  time \item sigma, a matrix (n x k) containing the prior
  standard deviations for the log-change in each lambda at
  each time \item a list phi, containing the numeric prior
  \code{df} and a numeric vector \code{scale} of length n
  }}

  \item{lambda0}{numeric vector (length k) of time 0 prior
  means for OD flows}

  \item{sigma0}{numeric vector (length k) of time 0 prior
  standard deviations for OD flows}

  \item{phi0}{numeric starting value for phi at time 0}

  \item{rho}{numeric fixed autoregressive parameter for
  dynamics on lambda; see reference for details}

  \item{tau}{numeric fixed power parameter for variance
  structure on truncated normal noise; see reference for
  details}

  \item{m}{integer number of particles to use}

  \item{verbose}{logical activates verbose diagnostic
  output}

  \item{Xdraws}{integer number of draws to perform for
  \code{xsample} RDA}

  \item{Xburnin}{integer number of burnin draws to discard
  for \code{xsample} proposals RDA in addition to baseline
  number of draws}

  \item{Movedraws}{integer number of iterations to run for
  each move step}

  \item{nThresh}{numeric effective number of independent
  particles below which redraw will be performed}

  \item{aggregate}{logical to activate aggregation of MCMC
  results; highly}

  \item{backward}{logical to activate reverse filtering
  (for smoothing}

  \item{tStart}{integer time index to begin iterations
  from}
}
\value{
list containing: \itemize{ \item xList \item lambdaList
\item phiList \item y \item rho \item prior \item n \item l
\item k \item A \item A_qr \item A1 \item A1_inv \item A2
\item nEff \item tStart \item backward \item aggregate }
}
\description{
Particle filtering with sample-resample-move algorithm for
multilevel state-space model of Blocker & Airoldi (2011).
This has log-normal autoregressive dynamics on OD
intensities, log-normal emission distributions, and
truncated normal observation densities.  This can return
full (all particles) output, but it is typically better to
aggregate results as you go to reduce memory consumption.
It can also run forward or backward filtering for
smoothing. These results are combined via a separate
function for smoothing; however, this procedure typically
performs poorly due to differences between the
distributions of particles from forward and reverse
filtering.
}
\references{
A.W. Blocker and E.M. Airoldi. Deconvolution of mixing time
series on a graph. Proceedings of the Twenty-Seventh
Conference Annual Conference on Uncertainty in Artificial
Intelligence (UAI-11) 51-60, 2011.
}
\seealso{
Other bayesianDynamicModel: \code{\link{buildPrior}};
\code{\link{move_step}}
}


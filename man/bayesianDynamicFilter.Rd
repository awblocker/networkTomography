\name{bayesianDynamicFilter}
\alias{bayesianDynamicFilter}
\title{Function for inference with multilevel state-space model}
\usage{
  bayesianDynamicFilter(Y, A, prior, lambda0, sigma0, phi0,
    rho = 0.1, tau = 2, m = 1000, verbose = FALSE,
    Xdraws = 5 * m, Xburnin = m, Movedraws = 10,
    nThresh = 10, aggregate = FALSE, backward = FALSE,
    tStart = 1)
}
\arguments{
  \item{Y}{}

  \item{A}{}

  \item{prior}{}

  \item{lambda0}{}

  \item{sigma0}{}

  \item{phi0}{}

  \item{rho}{}

  \item{tau}{}

  \item{m}{}

  \item{verbose}{}

  \item{Xdraws}{}

  \item{Xburnin}{}

  \item{Movedraws}{}

  \item{nThresh}{}

  \item{aggregate}{}

  \item{backward}{}

  \item{tStart}{}
}
\value{
  list containing: \itemize{ \item xList \item lambdaList
  \item phiList \item y \item rho \item prior \item n \item
  l \item k \item A \item A_qr \item A1 \item A1_inv \item
  A2 \item nEff \item tStart \item backward \item aggregate
  }
}
\description{
  (log-normal autoregressive dynamics, truncated normal
  observation densities) Can return full (all particles)
  output Can run forward or backward filtering; combination
  via separate function for smoothing
}
\references{
  A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
  time series on a graph. Proceedings of the Twenty-Seventh
  Conference Annual Conference on Uncertainty in Artificial
  Intelligence (UAI-11) 51-60, 2011.
}
\seealso{
  Other bayesianDynamicModel: \code{\link{move_step}}
}


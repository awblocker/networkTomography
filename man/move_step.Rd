\name{move_step}
\alias{move_step}
\title{Move step of sample-resample-move algorithm for multilevel state-space model}
\usage{
  move_step(y, X, tme, lambda, phi, lambdatm1, phitm1,
    prior, A, A1_inv, A2, rho, tau, m = ncol(X),
    l = nrow(A1_inv), k = length(lambda), ndraws = 10,
    minAccepts = 0, verbose = FALSE)
}
\arguments{
  \item{y}{}

  \item{X}{}

  \item{tme}{}

  \item{lambda}{}

  \item{phi}{}

  \item{lambdatm1}{}

  \item{phitm1}{}

  \item{prior}{}

  \item{A}{}

  \item{A1_inv}{}

  \item{A2}{}

  \item{rho}{}

  \item{tau}{}

  \item{m}{}

  \item{l}{}

  \item{k}{}

  \item{ndraws}{}

  \item{minAccepts}{}

  \item{verbose}{}
}
\value{
  list containing updated values of X, lambda, and phi
}
\description{
  Move step of sample-resample-move algorithm for
  multilevel state-space model
}
\references{
  A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
  time series on a graph. Proceedings of the Twenty-Seventh
  Conference Annual Conference on Uncertainty in Artificial
  Intelligence (UAI-11) 51-60, 2011.
}
\seealso{
  Other bayesianDynamicModel:
  \code{\link{bayesianDynamicFilter}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


\name{twMCMC}
\alias{twMCMC}
\title{Function to run MCMC sampling for model of Tebaldi & West (1998)}
\usage{
  twMCMC(Y, A, prior, ndraws = 120000, burnin = 20000,
    verbose = 0)
}
\arguments{
  \item{Y}{numeric vector of observed link loads at a
  single time (length k)}

  \item{A}{routing matrix of dimension (k x n); needs to be
  full row rank}

  \item{prior}{parameters for conjugate gamma prior
  (convolution and rate)}

  \item{ndraws}{integer number of draws for sampler to
  produce (excluding burn-in)}

  \item{burnin}{integer number of additional draws to
  discard as burnin}

  \item{verbose}{integer level of verbosity; levels > 1
  have no effect currently}
}
\value{
  list consisting of matrix of draws for X \code{XDraws},
  matrix of draws for X \code{lambdaDraws}, and vector of
  acceptances per OD flow \code{accepts}
}
\description{
  Runs MCMC sampling for the gamma-Poisson model presented
  in Tebaldi & West (1998). The algorithm used is a
  modification of that presented in the original paper. It
  uses a joint proposal for (x_k, lambda_k) to greatly
  accelerate convergence.
}
\examples{
cat("Example to come -- this is meant to annoy you into writing it",
    file=stderr())
}
\keyword{models}
\keyword{multivariate}


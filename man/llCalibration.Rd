\name{llCalibration}
\alias{llCalibration}
\title{Evaluate marginal log-likelihood for calibration SSM}
\usage{
llCalibration(theta, Ft, yt, Zt, Rt, k = ncol(Ft), tau = 2,
  initScale = 1/(1 - diag(Ft)^2), nugget = sqrt(.Machine$double.eps))
}
\arguments{
  \item{theta}{numeric vector (length k+1) of parameters.
  theta[-1] = log(lambda), and theta[1] = log(phi)}

  \item{Ft}{evolution matrix (k x k) for OD flows; include
  fixed}

  \item{yt}{matrix (k x n) of observed link loads, one
  observation per column}

  \item{Zt}{observation matrix for system; should be
  routing matrix A}

  \item{Rt}{covariance matrix for observation equation;
  typically small and fixed}

  \item{k}{integer number of OD flows to infer}

  \item{tau}{numeric power parameter for mean-variance
  relationship}

  \item{initScale}{numeric inflation factor for time-zero
  state covariance; defaults to steady-state variance
  setting}

  \item{nugget}{small positive value to add to diagonal of
  state evolution covariance matrix to ensure numerical
  stability}
}
\value{
numeric marginal log-likelihood obtained via Kalman
smoothing
}
\description{
Evaluates marginal log-likelihood for calibration SSM of
Blocker & Airoldi (2011) using Kalman filtering. This is
very fast and numerically stable, using the univariate
Kalman filtering and smoothing functions of \code{KFAS}
with Fortran implementations.
}
\references{
A.W. Blocker and E.M. Airoldi. Deconvolution of mixing time
series on a graph. Proceedings of the Twenty-Seventh
Conference Annual Conference on Uncertainty in Artificial
Intelligence (UAI-11) 51-60, 2011.
}
\seealso{
Other calibrationModel: \code{\link{calibration_ssm}};
\code{\link{mle_filter}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


\name{mle_filter}
\alias{mle_filter}
\title{Filtering & smoothing at MLE for calibration SSM}
\usage{
mle_filter(mle, Ft, yt, Zt, Rt, k = ncol(Ft), tau = 2, initScale = 1/(1 -
  diag(Ft)^2), nugget = sqrt(.Machine$double.eps))
}
\arguments{
  \item{mle}{numeric vector (length k+1) of parameters.
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

list containing result of Kalman smoothing; see
\code{\link{SSModel}} and \code{\link{KFS}} for details
}
\description{
Run Kalman filtering and smoothing at calculated MLE for
parameters of calibration SSM. This is used to obtain point
and covariance estimates for the actual OD flows X
following estimation of other parameters.
}
\references{
A.W. Blocker and E.M. Airoldi. Deconvolution of mixing time
series on a graph. Proceedings of the Twenty-Seventh
Conference Annual Conference on Uncertainty in Artificial
Intelligence (UAI-11) 51-60, 2011.
}
\seealso{
Other calibrationModel: \code{\link{calibration_ssm}};
\code{\link{llCalibration}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


\name{calibration_ssm}
\alias{calibration_ssm}
\title{Estimation for the linear SSM calibration model of Blocker & Airoldi (2011)}
\usage{
  calibration_ssm(tme, y, A, Ft, Rt, lambda0, phihat0,
    tau = 2, w = 11, initScale = 1/(1 - diag(Ft)^2),
    nugget = sqrt(.Machine$double.eps), verbose = FALSE,
    logTrans = TRUE, method = "L-BFGS-B",
    optimArgs = list())
}
\arguments{
  \item{tme}{integer time at which to center moving window
  for estimation}

  \item{y}{matrix (n x m) of observed link loads from all
  times (not just the window used for estimation; one
  observation per row}

  \item{A}{routing matrix (m x k) for network; should be
  full row rank}

  \item{Ft}{matrix (k x k) containing fixed autoregressive
  parameters for state evolution equation; upper-left block
  of overall matrix for expanded state}

  \item{Rt}{covariance matrix for observation equation;
  typically small and fixed}

  \item{lambda0}{matrix (n x k) of initial estimates for
  lambda (e.g. obtained via IPFP)}

  \item{phihat0}{numeric vector (length n) of initial
  estimates for phi}

  \item{tau}{numeric power parameter for mean-variance
  relationship}

  \item{w}{number of observations to use for rolling-window
  estimation; handles boundary cases cleanly}

  \item{initScale}{numeric inflation factor for time-zero
  state covariance; defaults to steady-state variance
  setting}

  \item{nugget}{small positive value to add to diagonal of
  state evolution covariance matrix to ensure numerical
  stability}

  \item{verbose}{logical to select verbose output from
  algorithm}

  \item{logTrans}{logical whether to log-transform
  parameters for optimization.  If FALSE, sets method to
  "L-BFGS-B".}

  \item{method}{optimization method to use (in optim
  calls)}

  \item{optimArgs}{list of arguments to append to control
  argument for optim.  Can include all arguments except for
  fnscale, which is automatically set}
}
\value{
  list containing \code{lambdahat}, a numeric vector
  (length k) containing the MLE for lambda; \code{phihat},
  the MLE for phi; \code{xhat}, the smoothed estimates of
  the OD flows for the window used as a k x w matrix; and
  \code{varhat}, a k x w matrix containing the diagonal of
  the estimated covariance for each OD flow in the window
}
\description{
  Maximum likelihood estimation of the parameters of the
  calibration model from Blocker & Airoldi (2011) via
  direct numerical maximization of the marginal
  log-likelihood. This relies upon efficient Kalman
  smoothing to evaluate the marginal likelihood, which is
  provided here by the \code{KFAS} package.
}
\examples{
data(bell.labs)

lambda0 <- matrix(1, nrow(bell.labs$Y), ncol(bell.labs$A))
lambda0[100,] <- ipfp(y=bell.labs$Y[100,], A=bell.labs$A,
                      x0=rep(1, ncol(bell.labs$A)))
phihat0 <- rep(1, nrow(bell.labs$Y))
Ft <- 0.5 * diag_mat(rep(1, ncol(bell.labs$A)))
Rt <- 0.01 * diag_mat(rep(1, nrow(bell.labs$A)))

# Not run
#fit.calibration <- calibration_ssm(tme=100, y=bell.labs$Y, A=bell.labs$A,
#                                   Ft=Ft, Rt=Rt, lambda0=lambda0,
#                                   phihat0=phihat0, w=23)
}
\references{
  A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
  time series on a graph. Proceedings of the Twenty-Seventh
  Conference Annual Conference on Uncertainty in Artificial
  Intelligence (UAI-11) 51-60, 2011.
}
\seealso{
  Other calibrationModel: \code{\link{llCalibration}},
  \code{\link{mle_filter}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


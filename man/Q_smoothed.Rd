\name{Q_smoothed}
\alias{Q_smoothed}
\title{Q function for smoothed EM algorithm of Cao et al. (2000)}
\usage{
Q_smoothed(logtheta, c, M, rdiag, eta0, sigma0, V, eps.lambda, eps.phi)
}
\arguments{
  \item{logtheta}{numeric vector (length k+1) of
  log(lambda) (1:k) and log(phi) (last entry)}

  \item{c}{power parameter in model of Cao et al. (2000)}

  \item{M}{matrix (n x k) of conditional expectations for
  OD flows, one time per row}

  \item{rdiag}{numeric vector (length k) containing
  diagonal of conditional covariance matrix R}

  \item{eta0}{numeric vector (length k+1) containing value
  for log(c(lambda, phi)) from previous time (or initial
  value)}

  \item{sigma0}{covariance matrix (k+1 x k+1) of
  log(c(lambda, phi)) from previous time (or initial
  value)}

  \item{V}{evolution covariance matrix (k+1 x k+1) for
  log(c(lambda, phi)) (random walk)}

  \item{eps.lambda}{numeric small positive value to add to
  lambda for numerical stability; typically 0}

  \item{eps.phi}{numeric small positive value to add to phi
  for numerical stability; typically 0}
}
\value{
numeric value of Q function; not vectorized in any way
}
\description{
Computes the Q function (expected log-likelihood) for the
EM algorithm of Cao et al. (2000) for their smoothed model.
}
\references{
J. Cao, D. Davis, S. Van Der Viel, and B. Yu. Time-varying
network tomography: router link data. Journal of the
American Statistical Association, 95:1063-75, 2000.
}
\seealso{
Other CaoEtAl: \code{\link{Q_iid}}; \code{\link{R_estep}};
\code{\link{grad_iid}}; \code{\link{grad_smoothed}};
\code{\link{locally_iid_EM}}; \code{\link{m_estep}};
\code{\link{phi_init}}; \code{\link{smoothed_EM}}
}


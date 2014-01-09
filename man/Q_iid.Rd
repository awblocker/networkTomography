\name{Q_iid}
\alias{Q_iid}
\title{Q function for locally IID EM algorithm of Cao et al. (2000)}
\usage{
Q_iid(logtheta, c, M, rdiag, epsilon)
}
\arguments{
  \item{logtheta}{numeric vector (length k+1) of
  log(lambda) (1:k) and log(phi) (last entry)}

  \item{c}{power parameter in model of Cao et al. (2000)}

  \item{M}{matrix (n x k) of conditional expectations for
  OD flows, one time per row}

  \item{rdiag}{numeric vector (length k) containing
  diagonal of conditional covariance matrix R}

  \item{epsilon}{numeric nugget to add to diagonal of
  covariance for numerical stability}
}
\value{
numeric value of Q function; not vectorized in any way
}
\description{
Computes the Q function (expected log-likelihood) for the
EM algorithm of Cao et al. (2000) for their locally IID
model.
}
\references{
J. Cao, D. Davis, S. Van Der Viel, and B. Yu. Time-varying
network tomography: router link data. Journal of the
American Statistical Association, 95:1063-75, 2000.
}
\seealso{
Other CaoEtAl: \code{\link{Q_smoothed}};
\code{\link{R_estep}}; \code{\link{grad_iid}};
\code{\link{grad_smoothed}}; \code{\link{locally_iid_EM}};
\code{\link{m_estep}}; \code{\link{phi_init}};
\code{\link{smoothed_EM}}
}


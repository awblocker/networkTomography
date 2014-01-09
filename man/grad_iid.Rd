\name{grad_iid}
\alias{grad_iid}
\title{Compute analytic gradient of Q-function for locally IID EM algorithm of Cao
et al. (2000)}
\usage{
grad_iid(logtheta, c, M, rdiag, epsilon)
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
numeric vector of same length as logtheta containing
calculated gradient
}
\description{
Computes gradient of Q-function with respect to
log(c(lambda,phi)) for EM algorithm from Cao et al. (2000)
for their locally IID model.
}
\references{
J. Cao, D. Davis, S. Van Der Viel, and B. Yu. Time-varying
network tomography: router link data. Journal of the
American Statistical Association, 95:1063-75, 2000.
}
\seealso{
Other CaoEtAl: \code{\link{Q_iid}};
\code{\link{Q_smoothed}}; \code{\link{R_estep}};
\code{\link{grad_smoothed}}; \code{\link{locally_iid_EM}};
\code{\link{m_estep}}; \code{\link{phi_init}};
\code{\link{smoothed_EM}}
}


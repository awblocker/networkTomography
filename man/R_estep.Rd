\name{R_estep}
\alias{R_estep}
\title{Compute conditional covariance matrix for EM algorithms of Cao et al. (2000)}
\usage{
  R_estep(lambda, phi, A, c, epsilon)
}
\arguments{
  \item{lambda}{numeric vector (length k) of mean OD flows
  from last M-step}

  \item{phi}{numeric scalar scale for covariance matrix of
  xt}

  \item{A}{routing matrix (m x k) for network being
  analyzed}

  \item{c}{power parameter in model of Cao et al. (2000)}

  \item{epsilon}{numeric nugget to add to diagonal of
  covariance for numerical stability}
}
\value{
  conditional covariance matrix (k x k) of OD flows given
  parameters
}
\description{
  Computes conditional covariance of OD flows for E-step of
  EM algorithm from Cao et al. (2000) for their locally IID
  model.
}
\references{
  J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
  Time-varying network tomography: router link data.
  Journal of the American Statistical Association,
  95:1063-75, 2000.
}
\seealso{
  Other CaoEtAl: \code{\link{grad_iid}},
  \code{\link{grad_smoothed}},
  \code{\link{locally_iid_EM}}, \code{\link{m_estep}},
  \code{\link{phi_init}}, \code{\link{Q_iid}},
  \code{\link{Q_smoothed}}, \code{\link{smoothed_EM}}
}


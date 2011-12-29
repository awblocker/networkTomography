\name{phi_init}
\alias{phi_init}
\title{Simple initialization for phi in model of Cao et al. (2000)}
\usage{
  phi_init(Y, A, lambda0, c)
}
\arguments{
  \item{Y}{matrix (n x k) of observed link loads over time}

  \item{A}{routing matrix (m x k)}

  \item{lambda0}{numeric vector (length k) of initial
  guesses for lambda}

  \item{c}{power parameter in model of Cao et al. (2000)}
}
\value{
  numeric starting value for phi
}
\description{
  Uses a crude estimator to get a starting point for phi in
  the model of Cao et al. (2000).
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
  \code{\link{Q_iid}}, \code{\link{Q_smoothed}},
  \code{\link{R_estep}}, \code{\link{smoothed_EM}}
}


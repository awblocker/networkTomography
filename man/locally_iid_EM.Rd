\name{locally_iid_EM}
\alias{locally_iid_EM}
\title{Run EM algorithm to obtain MLE for locally IID model of Cao et al. (2000)}
\usage{
  locally_iid_EM(Y, A, lambda0, phi0 = NULL, c = 2,
    maxiter = 1000, tol = 1e-06, epsilon = 0.01,
    method = "L-BFGS-B", checkActive = FALSE)
}
\arguments{
  \item{Y}{matrix (h x k) of observations in local window;
  columns correspond to OD flows, and rows are individual
  observations}

  \item{A}{routing matrix (m x k) for network being
  analyzed}

  \item{lambda0}{initial vector of values (length k) for
  lambda; \code{ipfp} is a good way to obtain this}

  \item{phi0}{initial value for covariance scale phi;
  initializes automatically using \code{phi_init} if NULL,
  but you can likely do better}

  \item{c}{power parameter in model of Cao et al. (2000)}

  \item{maxiter}{maximum number of EM iterations to run}

  \item{tol}{tolerance (in relative change in Q function
  value) for stopping EM iterations}

  \item{epsilon}{numeric nugget to add to diagonal of
  covariance for numerical stability}

  \item{method}{optimization method to use (in optim
  calls)}

  \item{checkActive}{logical check for deterministically
  known OD flows}
}
\value{
  list with 3 elements: \code{lambda}, the estimated value
  of lambda; \code{phi}, the estimated value of phi; and
  \code{iter}, the number of iterations run
}
\description{
  Runs EM algorithm to compute MLE for the locally IID
  model of Cao et al. (2000). Uses numerical optimization
  of Q-function for each M-step with analytic computation
  of its gradient.
}
\references{
  J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
  Time-varying network tomography: router link data.
  Journal of the American Statistical Association,
  95:1063-75, 2000.
}
\seealso{
  Other CaoEtAl: \code{\link{grad_iid}},
  \code{\link{grad_smoothed}}, \code{\link{m_estep}},
  \code{\link{phi_init}}, \code{\link{Q_iid}},
  \code{\link{Q_smoothed}}, \code{\link{R_estep}},
  \code{\link{smoothed_EM}}
}


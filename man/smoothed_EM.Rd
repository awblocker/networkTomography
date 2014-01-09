\name{smoothed_EM}
\alias{smoothed_EM}
\title{Run EM algorithm to obtain MLE (single time) for smoothed model of Cao et al.
(2000)}
\usage{
smoothed_EM(Y, A, eta0, sigma0, V, c = 2, maxiter = 1000, tol = 1e-06,
  eps.lambda = 0, eps.phi = 0, method = "L-BFGS-B")
}
\arguments{
  \item{Y}{matrix (h x k) of observations in local window;
  columns correspond to OD flows, and rows are individual
  observations}

  \item{A}{routing matrix (m x k) for network being
  analyzed}

  \item{eta0}{numeric vector (length k+1) containing value
  for log(c(lambda, phi)) from previous time (or initial
  value)}

  \item{sigma0}{covariance matrix (k+1 x k+1) of
  log(c(lambda, phi)) from previous time (or initial
  value)}

  \item{V}{evolution covariance matrix (k+1 x k+1) for
  log(c(lambda, phi)) (random walk)}

  \item{c}{power parameter in model of Cao et al. (2000)}

  \item{maxiter}{maximum number of EM iterations to run}

  \item{tol}{tolerance (in relative change in Q function
  value) for stopping EM iterations}

  \item{eps.lambda}{numeric small positive value to add to
  lambda for numerical stability; typically 0}

  \item{eps.phi}{numeric small positive value to add to phi
  for numerical stability; typically 0}

  \item{method}{optimization method to use (in optim
  calls)}
}
\value{
list with 5 elements: \code{lambda}, the estimated value of
lambda; \code{phi}, the estimated value of phi;
\code{iter}, the number of iterations run; \code{etat},
log(c(lambda, phi)); and sigmat, the inverse of the Q
functions Hessian at its mode
}
\description{
Runs EM algorithm to compute MLE for the smoothed model of
Cao et al. (2000). Uses numerical optimization of
Q-function for each M-step with analytic computation of its
gradient. This performs estimation for a single time point
using output from the previous one.
}
\references{
J. Cao, D. Davis, S. Van Der Viel, and B. Yu. Time-varying
network tomography: router link data. Journal of the
American Statistical Association, 95:1063-75, 2000.
}
\seealso{
Other CaoEtAl: \code{\link{Q_iid}};
\code{\link{Q_smoothed}}; \code{\link{R_estep}};
\code{\link{grad_iid}}; \code{\link{grad_smoothed}};
\code{\link{locally_iid_EM}}; \code{\link{m_estep}};
\code{\link{phi_init}}
}


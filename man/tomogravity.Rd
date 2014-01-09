\name{tomogravity}
\alias{tomogravity}
\title{Run tomogravity estimation on complete time series of aggregate flows}
\usage{
tomogravity(Y, A, lambda, lower = 0, normalize = FALSE,
  .progress = "none", control = list())
}
\arguments{
  \item{Y}{n x m matrix contain one vector of observed
  aggregate flows per row.  This should include all
  observed aggegrate flows with none removed due to
  redundancy.}

  \item{A}{m x k routing matrix. This need not be of full
  row rank and must include all source and destination
  flows.}

  \item{lambda}{Regularization parameter for mutual
  information prior. Note that this is scaled by the
  squared total traffic in the objective function before
  scaling the mututal information prior.}

  \item{lower}{Component-wise lower bound for xt in
  L-BFGS-B optimization.}

  \item{normalize}{If TRUE, xt and yt are scaled by N.
  Typically used in conjunction with calcN to normalize
  traffic to proportions, easing the tuning of lambda.}

  \item{.progress}{name of the progress bar to use, see
  \code{\link{create_progress_bar}} in plyr documentation}

  \item{control}{List of control information for optim.}
}
\value{
A list containing three elements: \itemize{ \item
resultList, a list containing the output from running
\code{\link{tomogravity.fit}} on each timepoint \item
changeFromInit, a vector of length n containing the
relative L_1 change between the initial (IPFP)
point-to-point flow estimates and the final tomogravity
estimates \item Xhat, a n x k matrix containing a vector of
estimated point-to-point flows (for each time point) per
row }
}
\description{
The aggregate flows Y and their corresponding routing
matrix A must include all aggregate source and destination
flows.
}
\examples{
data(cmu)
estimate <- tomogravity(Y=cmu$Y.full[1, , drop=FALSE], A=cmu$A.full,
                        lambda=0.01, .progress='text')
}
\seealso{
Other tomogravity: \code{\link{tomogravity.fit}}
}
\keyword{models}
\keyword{multivariate}
\keyword{ts}


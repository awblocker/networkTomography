
#' Evaluate marginal log-likelihood for calibration SSM
#'
#' Evaluates marginal log-likelihood for calibration SSM of Blocker & Airoldi
#' (2011) using Kalman smoothing. This is very fast and numerically-stable,
#' using the univariate Kalman filtering and smoothing functions of \code{KFAS}
#' with Fortran implementations.
#'
#' @param theta numeric vector (length k+1) of parameters. theta[-1] =
#'      log(lambda), and theta[1] = log(phi)
#' @param Ft evolution matrix (2k x 2k) for OD flows; include fixed
#       autoregressive parameters
#' @param qind integer vector of indices to update lambda portion of Ft matrix;
#'      typically a subset of \code{diag_ind(2*k)}
#' @param yt matrix (k x n) of observed link loads, one observation per column
#' @param Zt observation matrix for system; should be padded version of routing
#'      matrix A
#' @param R covariance matrix for observation equation; typically small and
#'      fixed
#' @param k integer number of OD flows to infer
#' @param tau numeric power parameter for mean-variance relationship
#' @param scale numeric inflation factor for time-zero state covariance
#' @param nugget small positive value to add to diagonal of state evolution
#'      covariance matrix to ensure numerical stability
#' @return numeric marginal log-likelihood obtained via Kalman smoothing
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family calibrationModel
llCalibration <- function(theta, Ft, qind, yt, Zt, R,
                         k=ncol(Ft)/2, tau=2, scale=10,
                         nugget=sqrt(.Machine$double.eps)) {
    # Parse parameters
    lambda <- exp(theta[-1])
    phi <- exp(theta[1])

    P1 <- diag_mat(scale*c(phi*lambda^tau + nugget,rep(0,k)))
    a1 <- c(lambda, rep(1,k))

    # Setup matrices
    Ft[qind] <- lambda
    V <- diag_mat(c(phi*lambda^tau + nugget, rep(0,k)))

    # Run Kalman filter & smoother
    f.out <- kf(yt=yt, Zt=Zt, Tt=Ft, Rt=diag(ncol(Ft)), Ht=R, Qt=V,
                a1=a1, P1=P1)
    f.out <- ks(f.out)

    # Return log-likelihood
    return(f.out$lik)
}

#' Filtering & smoothing at MLE for calibration SSM
#'
#' Run Kalman filtering and smoothing at calculated MLE for parameters of
#' calibration SSM. This is used to obtain point and covariance estimates for
#' the actual OD flows X following estimation.
#'
#' @param mle numeric vector (length k+1) of parameters. theta[-1] =
#'      log(lambda), and theta[1] = log(phi)
#' @param Ft evolution matrix (2k x 2k) for OD flows; include fixed
#       autoregressive parameters
#' @param qind integer vector of indices to update lambda portion of Ft matrix;
#'      typically a subset of \code{diag_ind(2*k)}
#' @param yt matrix (k x n) of observed link loads, one observation per column
#' @param Zt observation matrix for system; should be padded version of routing
#'      matrix A
#' @param R covariance matrix for observation equation; typically small and
#'      fixed
#' @param k integer number of OD flows to infer
#' @param tau numeric power parameter for mean-variance relationship
#' @param scale numeric inflation factor for time-zero state covariance
#' @param nugget small positive value to add to diagonal of state evolution
#'      covariance matrix to ensure numerical stability
#' @return numeric marginal log-likelihood obtained via Kalman smoothing
#' @return list containing result of Kalman smoothing; see \code{\link{ks}} for 
#'      details
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family calibrationModel
mle_filter <- function(mle, Ft, qind, yt, Zt, R,
                       k=ncol(Ft)/2, tau=2, scale=10,
                       nugget=sqrt(.Machine$double.eps)) {
    # Parse parameters
    lambda <- exp(mle$par[-1])
    phi <- exp(mle$par[1])

    P1 <- diag_mat(scale*c(phi*lambda^tau + nugget,rep(0,k)))
    a1 <- c(lambda, rep(1,k))

    # Setup matrices
    Ft[qind] <- lambda
    V <- diag_mat(c(phi*lambda^tau + nugget, rep(0,k)))

    # Run Kalman filter & smoother
    f.out <- kf(yt=yt, Zt=Zt, Tt=Ft, Rt=diag(ncol(Ft)), Ht=R, Qt=V,
                a1=a1, P1=P1)
    f.out <- ks(f.out)
}

#' Estimation for the linear SSM calibration model of Blocker & Airoldi (2011)
#'
#' Maximum likelihood estimation of the parameters of the calibration model from
#' Blocker & Airoldi (2011) via direct numerical maximization of the marginal
#' log-likelihood. This relies upon efficient Kalman smoothing to evaluate the
#' marginal likelihood, which is provided here by the \code{KFAS} package.
#'
#' @param tme integer time at which to center moving window for estimation
#' @param y matrix (n x m) of observed link loads from all times (not just the
#'      window used for estimation; one observation per row
#' @param A routing matrix (m x k) for network; should be full row rank
#' @param F matrix (k x k) containing fixed autoregressive parameters for state
#'      evolution equation; upper-left block of overall matrix for expanded
#'      state
#' @param R covariance matrix for observation equation; typically small and
#'      fixed
#' @param xhat0 matrix (n x k) of initial estimates for OD flows X (e.g.
#'      obtained via IPFP)
#' @param phihat0 numeric vector (length n) of initial estimates for phi
#' @param tau numeric power parameter for mean-variance relationship
#' @param w number of observations to use for rolling-window estimation; handles
#'      boundary cases cleanly
#' @param maxiter maximum number of iterations to use in numerical optimization
#'      of log-likelihood
#' @param tol tolerance to use for numerical optimization
#' @param scale numeric inflation factor for time-zero state covariance
#' @param nugget small positive value to add to diagonal of state evolution
#'      covariance matrix to ensure numerical stability
#' @param verbose logical to select verbose output from algorithm
#' @param method optimization method to use (in optim calls)
#' @return list containing \code{lambdahat}, a numeric vector (length k)
#'      containing the MLE for lambda; \code{phihat}, the MLE for phi;
#'      \code{xhat}, the smoothed estimates of the OD flows for the window used
#'      as a k x w matrix; and \code{varhat}, a k x w matrix containing the 
#'      diagonal of the estimated covariance for each OD flow in the window
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family calibrationModel
calibration_ssm <- function(tme, y, A, F, R, xhat0, phihat0,
                            tau=2, w=11, maxiter=1e4, tol=1e-9, scale=10,
                            nugget=sqrt(.Machine$double.eps), verbose=FALSE,
                            method='Nelder-Mead') {
    # Calculate dimensions
    k <- ncol(A)
    l <- ncol(y)

    # Calculate window parameters
    h <- floor(w/2)

    # Calculate length of window and index of tme within window
    if ( (tme-h>0) & (tme+h<nrow(y)) ) {
        n <- w
        t_ind <- h+1
    } else {
        # Handle border cases
        if (tme-h<=0) {
            n <- tme+h
            t_ind <- tme
        } else {
            n <- nrow(y)-tme+h+1
            t_ind <- h+1
        }
    }

    # Print tme and window length if verbose
    if (verbose) {
        cat(sprintf('tme = %d; n = %d\n', tme, n), file=stderr())
    }

    # Setup array for F
    Ft <- matrix(0, 2*k, 2*k)
    dind <- diag_ind(2*k)
    # Diagonal block
    Ft[dind[seq(k+1,2*k)]] <- 1
    # F block
    Ft[1:k, 1:k] <- F

    # Indices for Q block
    qind <- 2L*k*k + 1 + 0L:(k-1L)*(2L*k+1L)

    # Setup Zt
    Zt <- matrix(0, l, 2*k)
    Zt[,1:k] <- A

    # Setup data structures for parameters of interest
    x <- xhat0[max(tme-h,1):min(tme+h,nrow(y)),]

    # Setup data
    yt <- t( y[max(1,tme-h):min(tme+h,nrow(y)),] )
    lambda <- colMeans(xhat0[max(1,tme-h):min(tme+h,nrow(y)),])
    phi <- mean(phihat0[max(1,tme-h):min(tme+h,nrow(y))])

    # Setup starting values
    # a1 <- c(xhat0[max(1,tme-h-1),], rep(1,k))

    # Print starting values if verbose
    if (verbose) {
        cat('Starting value for lambda:\n', file=stderr())
        print(lambda)
        cat('Starting value for phi:\n', file=stderr())
        print(phi)
    }

    # Run numerical optimization using prediction error formulation
    # of Airoldi 2003 (section 3.1.2); with a Fortran implementation
    # of the univariate Kalman filter and smoother of Koopman & Durbin (2000,
    # 2003), this is far more efficient than the EM iterations (quadratic vs.
    # linear convergence, plus much less memory usage in the smoothing stage).
    mle <- optim(par=log(c(phi, lambda)), fn=llCalibration,
                 Ft=Ft, qind=qind, yt=yt, Zt=Zt, R=R, tau=tau, scale=scale,
                 nugget=nugget,
                 method=method,
                 control=list(fnscale=-1))

    # Print optim diagnostics if verbose
    if (verbose) {
        cat(sprintf('Convergence code: %d\n', mle$convergence), file=stderr())
        cat('Function evaluations:\n', file=stderr())
        print(mle$counts)
    }

    # Obtain Kalman filter output at MLE
    f.out <- mle_filter(mle=mle, Ft=Ft, qind=qind, yt=yt, Zt=Zt, R=R, tau=tau,
                        scale=scale, nugget=nugget)
    varhat <- apply(f.out$Pt, 3, diag)

    # Return results
    return(list(lambdahat=exp(mle$par[-1]), phihat=exp(mle$par[1]),
                xhat=f.out$ahat[1:k,t_ind], varhat=varhat[1:k,t_ind]))
}


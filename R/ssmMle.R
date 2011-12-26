
# Function for negative log-likelihood evaluation for SSM calibration
obj <- function(theta, Ft, qind, yt, Zt, R,
    k=ncol(Ft)/2, tau=2, scale=10, nugget=sqrt(.Machine$double.eps)) {
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

    # Return negative log-likelihood
    return(-f.out$lik)
}

# Function for filtering & smoothing at MLE of calibration SSM
mle_filter <- function(mle, Ft, qind, yt, Zt, R,
    k=ncol(Ft)/2, tau=2, scale=10, nugget=sqrt(.Machine$double.eps)) {
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

# Function for estimation of the linear SSM calibration model
# of Blocker & Airoldi (2011)
calibration_ssm <- function(tme, y, A, F, R, xhat0, phihat0,
    tau=2, w=11, maxiter=1e4, tol=1e-9, scale=10,
    nugget=sqrt(.Machine$double.eps), verbose=FALSE) {
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
        cat(sprintf('tme = %d; n = %d\n', tme, n))
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
        cat('Starting value for lambda:\n')
        print(lambda)
        cat('Starting value for phi:\n')
        print(phi)
    }

    # Run numerical optimization using prediction error formulation
    # of Airoldi 2003 (section 3.1.2); with a Fortran implementation
    # of the univariate Kalman filter and smoother of Koopman & Durbin (2000,
    # 2003), this is far more efficient than the EM iterations (quadratic vs.
    # linear convergence, plus much less memory usage in the smoothing stage).
    mle <- optim( log(c(phi, lambda)), obj,
        Ft=Ft, qind=qind, yt=yt, Zt=Zt, R=R, tau=tau, scale=scale,
        nugget=nugget,
        # method='BFGS')
        method='Nelder-Mead')

    # Print optim diagnostics if verbose
    if (verbose) {
        cat(sprintf('Convergence code: %d\n', mle$convergence))
        cat('Function evaluations:\n')
        print(mle$counts)
    }
    
    # Obtain Kalman filter output at MLE
    f.out <- mle_filter( mle, Ft=Ft, qind=qind, yt=yt, Zt=Zt, R=R,
        tau=tau, scale=scale, nugget=nugget)
    varhat <- apply(f.out$Pt, 3, diag)

    # Return results
    return(list(lambdahat=exp(mle$par[-1]), phihat=exp(mle$par[1]),
        xhat=f.out$ahat[1:k,t_ind], varhat=varhat[1:k,t_ind]))
}


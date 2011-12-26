
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

# Function for move step of sample-resample-move algorithm for multilevel
# state-space model with log-Normal / truncated Normal distributions
# (standard parameterizations)
move_step <- function(y, X, tme, lambda, phi,
    lambdatm1, phitm1, prior,
    A, A1_inv, A2,
    rho, tau,
    m=ncol(X), l=nrow(A1_inv), k=length(lambda),
    ndraws=10, minAccepts=0,
    verbose=FALSE) {

    # Setup X1 and X2 matrices (want each particle in a column)
    X1 <- t(X[,seq(1,l)])
    X2 <- t(X[,seq(l+1,k)])
    X <- t(X)
    
    X1prop <- X1
    X2prop <- X2
    Xprop <- X

    # Setup lambda matrices
    lambda <- t(lambda)
    lambdatm1 <- t(lambdatm1)

    # Initialize information for acceptance rates
    xAccepts <- matrix(0,m,k-l)
    lambdaAccepts <- matrix(0,m,k)
    phiAccepts <- rep(0,m)

    # Setup matrices for mean, std dev, and truncation adjustment for truncated
    # normal
    meanMat <- exp(lambda)
    sdMat <- sqrt( t( t(meanMat)^tau * (exp(phi^2) - 1) ) )
    # sdMat <- t( t(meanMat) * sqrt(exp(phi^2) - 1) )
    truncMat <- matrix(
        pnorm(0, meanMat, sdMat, lower.tail=FALSE, log.p=TRUE),
        k, m)
    
    # Loop for MCMC iterations
    iter <- 0
    while (iter < ndraws || (min(xAccepts) < minAccepts)) {
        # Draw lambda, phi | X via Gibbs steps

        # MH step on lambda
        varPropMat <- 1 / (1/prior$sigma[tme,]^2 +
            t(matrix(1/phi^2, m, k)))
        # meanPropMat <- ( ( prior$mu[tme,]+rho*lambdatm1 ) / prior$sigma[tme,]^2
        #     + t( t(log(X))/phi^2 ) ) * varPropMat
        meanPropMat <- lambda
        
        lambdaProp <- matrix(rnorm(k*m, meanPropMat, sqrt(varPropMat)), k, m)
        
        varPropMatRev <- varPropMat
        # meanPropMatRev <- meanPropMat
        meanPropMatRev <- lambdaProp
        
        # Calculate LLR and LIR
        meanMatProp <- exp(lambdaProp)
        # sdMatProp <- t( t(meanMatProp) * sqrt(exp(phi^2) - 1) )
        sdMatProp <- sqrt( t( t(meanMatProp)^tau * (exp(phi^2) - 1) ) )
        truncMatProp <- matrix(
            pnorm(0, meanMatProp, sdMatProp, lower.tail=FALSE, log.p=TRUE),
            k, m)

        llr <- matrix(
            dnorm(lambdaProp, prior$mu[tme,] + rho*lambdatm1,
            prior$sigma[tme,], log=TRUE) +
            dnorm(X, meanMatProp, sdMatProp, log=TRUE) - truncMatProp,
            k, m)
        llr <- llr - matrix(
            dnorm(lambda, prior$mu[tme,] + rho*lambdatm1,
            prior$sigma[tme,], log=TRUE) +
            dnorm(X, meanMat, sdMat, log=TRUE) - truncMat,
            k, m)

        lir <- matrix(
            dnorm(lambdaProp, meanPropMat, sqrt(varPropMat), log=TRUE),
            k, m)
        lir <- lir - matrix(
            dnorm(lambda, meanPropMatRev, sqrt(varPropMatRev), log=TRUE),
            k, m)
        
        # Accept with correct probability
        logAcceptProb <- (llr - lir)
        logAcceptProb[is.na(logAcceptProb)] <- -Inf

        acceptMat <- ( log(runif(k*m)) < logAcceptProb )

        # Handle acceptances and store necessary information
        lambda[acceptMat] <- lambdaProp[acceptMat]
        lambdaAccepts <- lambdaAccepts + t(acceptMat)
        
        # Update means and variances for X
        meanMat <- exp(lambda)
        # meanMat <- t( exp(t(lambda) + phi^2/2) )
        # sdMat <- t( t(meanMat) * sqrt(exp(phi^2) - 1) )
        sdMat <- sqrt( t( t(meanMat)^tau * (exp(phi^2) - 1) ) )
        truncMat <- matrix(
            pnorm(0, meanMat, sdMat, lower.tail=FALSE, log.p=TRUE),
            k, m)

        ###

        # MH step on phi
        alphaProp <- (k/2 + prior$phi$df)
        betaProp <- colSums((log(X)-lambda)^2)/2 + prior$phi$scale[tme]
        phiProp <- 1 / rgamma(m, alphaProp, betaProp)
        phiProp <- sqrt(phiProp)
        
        # Calculate LLR & LIR
        meanMatProp <- exp(lambda)
        # sdMatProp <- t( t(meanMatProp) * sqrt(exp(phiProp^2) - 1) )
        sdMatProp <- sqrt( t( t(meanMatProp)^tau * (exp(phiProp^2) - 1) ) )
        truncMatProp <- matrix(
            pnorm(0, meanMatProp, sdMatProp, lower.tail=FALSE, log.p=TRUE),
            k, m)

        llr <- (
            dgamma(1/phiProp^2, prior$phi$df, prior$phi$scale[tme], log=TRUE) -
            2*log(phiProp^2) +
            colSums(matrix(
            dnorm(X, meanMatProp, sdMatProp, log=TRUE) - truncMatProp,
            k, m))
            )
        llr <- llr - (
            dgamma(1/phi^2, prior$phi$df, prior$phi$scale[tme], log=TRUE) -
            2*log(phi^2) +
            colSums(matrix(
            dnorm(X, meanMat, sdMat, log=TRUE) - truncMat,
            k, m))
            )

        lir <- ( dgamma(1/phiProp^2, alphaProp, betaProp, log=TRUE) -
            2*log(phiProp^2) )
        lir <- lir - ( dgamma(1/phi^2, alphaProp, betaProp, log=TRUE) -
            2*log(phi^2) )

        # Accept with correct probability
        logAcceptProb <- (llr - lir)
        logAcceptProb[is.na(logAcceptProb)] <- -Inf

        acceptVec <- ( log(runif(m)) < logAcceptProb )

        # Handle acceptances and store necessary information
        phi[acceptVec] <- phiProp[acceptVec]
        phiAccepts <- phiAccepts + acceptVec

        # Update means and variances for X
        meanMat <- exp(lambda)
        # meanMat <- t( exp(t(lambda) + phi^2/2) )
        sdMat <- sqrt( t( t(meanMat)^tau * (exp(phi^2) - 1) ) )
        # sdMat <- t( t(meanMat) * sqrt(exp(phi^2) - 1) )
        truncMat <- matrix(
            pnorm(0, meanMat, sdMat, lower.tail=FALSE, log.p=TRUE),
            k, m)

        ###

        # Draw X | lambda, phi

        # # Using coordinate direction algorithm of Smith (1984)
        # # Draw (uniform) random coordinate, then propose uniformly along
        # Using random direction algorithm of Smith (1984)
        # Draw (uniform) random direction, then propose uniformly along
        # feasible line in that direction
        # Symmetric proposal -> random walk MH
        
        for (j in 1:(k-l)) {
            # # Calculate limits of feasible region along sampled directions
            # j <- sample(k-l,1)
            # remainderMat <- A1_inv %*% (y - A2[,-j] %*% X2[-j,])
            # adjVec <- as.numeric(A1_inv %*% A2[,j])
            # maxVec <- apply(remainderMat/adjVec, 2, function(col)
            #     min(col[adjVec>0]))
            # minVec <- apply(remainderMat/adjVec, 2, function(col)
            #     max(col[adjVec<0]))
            # minVec <- pmax(0,minVec)

            # Draw random directions (uniform on unit sphere)
            # by normalizing MV normal RVs
            dMat <- matrix( rnorm((k-l)*m), k-l, m )
            dMat <- t( t(dMat) / sqrt(colSums(dMat*dMat)) ) 
            
            # Calculate relevant matrices for feasible region calculations
            remainderMat <- X1
            adjMat <- A1_inv %*% A2 %*% dMat
            interceptMat <- remainderMat / adjMat
            
            # Find vector of maxima for \delta X2 st X1 >= 0
            maxVec1 <- sapply(1:m, function(i)
                if (max(interceptMat[,i])>0) 
                    return(min(interceptMat[,i][adjMat[,i]>0]))
                else
                    return(0)
                )
            # Find vector of maxima for \delta X2 st X2 >= 0
            maxVec2 <- sapply(1:m, function(i)
                if (max(-X2[,i]/dMat[,i])>0) 
                    return(min( (-X2[,i]/dMat[,i])[-dMat[,i]>0] ))
                else
                    return(0)
                )
            # Merge maxima
            maxVec <- pmin(maxVec1, maxVec2)

            # Find vector of minima for \delta X2 st X1 >= 0
            minVec1 <- sapply(1:m, function(i)
                if (min(interceptMat[,i])<0)
                    return(max(interceptMat[,i][adjMat[,i]<0]))
                else
                    return(0)
                )
            # Find vector of minima for \delta X2 st X2 >= 0
            minVec2 <- sapply(1:m, function(i)
                if (min(-X2[,i]/dMat[,i])<0)
                    return(max( (-X2[,i]/dMat[,i])[-dMat[,i]<0] ))
                else
                    return(0)
                )
            # Merge minima
            minVec <- pmax(minVec1, minVec2)

            proposal <- runif(m, minVec, maxVec)
            # pLower <- pnorm(minVec, meanMat[l+j,], sdMat[l+j,])
            # pUpper <- pnorm(maxVec, meanMat[l+j,], sdMat[l+j,])
            # u <- runif(m, pLower, pUpper)
            # proposal <- qnorm(u, meanMat[l+j,], sdMat[l+j,])
            
            # X2prop <- X2
            # X2prop[j,] <- proposal
            # X1prop <- A1_inv %*% (y - A2 %*% X2prop)
            X2prop <- X2 + t( t(dMat) * proposal )
            X1prop <- X1 - t( t(adjMat) * proposal )
            #Xprop <- rbind(X1prop, X2prop)
            Xprop[1:l,] <- X1prop
            Xprop[(l+1):k,] <- X2prop

            llr <- colSums( matrix(
                dnorm(Xprop, meanMat, sdMat, log=TRUE) -
                dnorm(X, meanMat, sdMat, log=TRUE),
                k, m) )
            lir <- 0
            # lir <-(
            #     dnorm(proposal, meanMat[l+j,], sdMat[l+j,], log=TRUE) -
            #     dnorm(X2[j,], meanMat[l+j,], sdMat[l+j,], log=TRUE)
            #     )

            logAcceptProb <- rep(-Inf, m)
            validVec <- !is.na(llr-lir)
            logAcceptProb[validVec] <- (llr - lir)[validVec]

            acceptVec <- (log(runif(m)) < logAcceptProb)

            # X2[j,acceptVec] <- proposal[acceptVec]
            # X1 <- A1_inv %*% (y - A2 %*% X2)
            # X <- rbind(X1, X2)
            X[,acceptVec] <- Xprop[,acceptVec]

            xAccepts[,j] <- xAccepts[,j] + acceptVec
        }
        
        iter <- iter + 1
    }

    # If verbose, print acceptance rate summaries
    if (verbose) {
        print(iter)
        # print(colMeans(xAccepts))
        print(summary(rowSums(xAccepts)))
        print(colMeans(lambdaAccepts))
        print(mean(phiAccepts))
    }
    
    # Return updated results
    return(list( X=t(X), lambda=t(lambda), phi=phi ))
}

# Utility functions for reparameterized log-normal distribution
# Reparameterization is based on Airoldi & Faloutsos 2004:
# X ~ logN( lambda, phi; tau ) -> E(X) = lambda, Var(X) = lambda^tau * phi

# RNG
rlnorm_mv <- function(n, lambda, phi, tau=2) {
    sigma <- sqrt(log(1+phi*lambda^(tau-2)))
    mu <- log(lambda) - sigma^2/2
    return( rlnorm(n, mu, sigma) )
}

# Density
dlnorm_mv <- function(x, lambda, phi, tau=2, logp=FALSE) {
    sigma <- sqrt(log(1+phi*lambda^(tau-2)))
    mu <- log(lambda) - sigma^2/2
    return( dlnorm(x, mu, sigma, log=logp) )
}

# Function to aggregate results; defaults to mean, SD, limits, and
# given quantiles
#
# Needed to prevent out-of-memory issues
agg <- function(mat, q=c(0.05, 0.16, 0.5, 0.84, 0.95)) {
    # Convert to matrix if needed
    if (is.vector(mat)) {
        mat <- as.matrix(mat)
    }
    
    ans <- matrix(0, length(q)+4, ncol(mat))
    nc <- ncol(mat)
    
    # Dimnames
    dimnames(ans) <- list()
    dimnames(ans)[[1]] <- vector("character", nrow(ans))
    dimnames(ans)[[1]][1:4] <- c("mean", "sd", "min", "max")
    
    # Default aggregates
    ans[1,] <- colMeans(mat)
    ans[2,] <- sd(mat)
    ans[3:4,] <- sapply( 1:nc, function(j) range(mat[,j]) )
    
    # Quantiles
    if (length(q) > 0) {
        ans[5:nrow(ans),] <- sapply( 1:nc, function(j)
            quantile(mat[,j], q) )
        dimnames(ans)[[1]][5:nrow(ans)] <- sprintf("q%02d", q*100)
    }
    
    return(ans)
}

# Function for inference with multilevel state-space model
# (log-normal autoregressive dynamics, truncated normal observation densities)
# Can return full (all particles) output
# Can run forward or backward filtering; combination via seperate function for
# smoothing
bayesianDynamicFilter <- function (Y, A, prior,
    lambda0, sigma0, phi0, 
    rho=0.1, tau=2, m=1e3,
    verbose=FALSE, # save.X=FALSE,
    Xdraws=5*m, Xburnin=m, Movedraws=10,
    nThresh=10,
    aggregate=FALSE,
    backward=FALSE,
    tStart=1) {
    # Calculate dimensions
    k <- ncol(A)
    l <- nrow(A)
    n <- nrow(Y)

    # Decompose A matrix into full-rank and singular parts; retain pivot info
    A_qr <- qr(A)
    pivot <- A_qr$pivot
    A_pivot <- A[,pivot]
    A1 <- A_pivot[,seq(A_qr$rank)]
    A1_inv <- solve(A1)
    A2 <- A_pivot[,seq(A_qr$rank+1,k)]

    # If verbose, print pivot information
    if (verbose) {
        cat('Pivot information:\n')
        print(pivot)
        cat(sprintf('rank(A) = %d\n', A_qr$rank))
    }

    # Reconfigure prior based on pivot (only using components pertaining to
    # X2, as X1 is determined based on data)
    prior$mu <- prior$mu[,pivot]
    prior$sigma <- prior$sigma[,pivot]

    # Reconfigure initial values based on pivot
    lambda0 <- lambda0[pivot]
    sigma0 <- sigma0[pivot]

    # If running backwards filter, reverse relevant matrices
    # Assumes prior has been input reversed; only reversing y
    if (backward) {
        Y <- Y[n:1,]
    }

    # Setup data structure for estimated OD flows
    xList <- list()
    lambdaList <- list()
    phiList <- list()
    nEffVec <- rep(0, n)

    # Setup data structures for particle filter
    lambda <- matrix(NA, m, k)
    w <- numeric(m)
    
    # Initialize particles from prior at time 0
    lambda <- matrix( rnorm(k*m, lambda0, sigma0), m, k, byrow=TRUE)
    phi <- rep(phi0, m)
    
    # Iterate over times
    for (tme in tStart:n) {
        # If verbose, print time
        if (verbose) {
            cat(sprintf('t = %d\n', tme))
            startTime <- proc.time()
        }

        # Retain last set of particles for move step
        lambdatm1 <- lambda
        phitm1 <- phi

        # # Draw set of phis from ad-hoc proposal
        # phi_prop <- rlnorm(m, log(phitm1)-phiscale^2/2, phiscale)
        # Draw phi from prior
        phi_prop <- 1/rgamma(m, prior$phi$df, prior$phi$scale[tme])
        phi_prop <- sqrt(phi_prop)
        
        # Draw set of lambdas
        lambda_prop <- ( rho*lambdatm1 +
            matrix(rnorm(k*m, prior$mu[tme,], prior$sigma[tme,]),
            m, k, byrow=TRUE) )
        
#        # Initialize X draws via Betrand's HNF method
#        x0 <- bhaas$samplesol2(Y[tme,], A_pivot, 2, alpha=1)$sampleMat[1,]

        # Draw X's from truncated normal on feasible region via xsample, 
        # then thin
        propMean <- exp( rho*colMeans(lambdatm1) + prior$mu[tme,] )
            # + (prior$phi$scale[tme]/prior$phi$df)/2)
        # propSD <- sqrt(exp(prior$phi$scale[tme]/prior$phi$df)-1)*propMean
        propSD <- sqrt( (exp(prior$phi$scale[tme]/prior$phi$df)-1) *
            propMean^tau )
        
        # Check for zero link flows
        activeLink <- !(Y[tme,] == 0)
        # activeOD <- !pmin( colSums(A_pivot[!activeLink,]), 1 )
        odRanges <- xranges(E=A_pivot, F=Y[tme,], ispos=TRUE)
        activeOD <- (odRanges[,2]-odRanges[,1]>0)
        nActive <- sum(activeOD)
        
        if (verbose) {
            cat(sprintf("nActive = %d\n", nActive))
        }
        
        # Start with IPF beginning from propMean
        x0 <- ipf(Y[tme,], A_pivot, propMean, tol=1e-6, maxit=Xdraws*1e2,
            verbose=FALSE)
        
        # Now, lsei for refinement
        x0Active <- lsei(
            A=diag(nActive), B=x0[activeOD],
            E=A_pivot[activeLink,activeOD,drop=FALSE], F=Y[tme,activeLink],
            G=diag(nActive), H=rep(0,nActive),
            tol=sqrt(.Machine$double.eps),
            type=1)$X
        
#        # Start with mirror algorithm
#        # Wrong distribution (uniform, not normal) and slow, but
#        # gets out of corner that lsei found
#        x0Active <- xsample(
#            E=A_pivot[activeLink,activeOD], F=Y[tme,activeLink],
#            G=diag(nActive), H=rep(0,nActive),
#            iter=10,
#            x0=x0Active, type='mirror')$X
#        x0Active <- drop(tail(x0Active, 1))
#        
#        if (verbose) {
#            cat(sprintf("Mirror complete; starting RDA\n"))
#        }
        
        # Switch to RDA algorithm
        # Correct distribution (truncated normal) and faster
        X_prop_active <- xsample(
            A=diag(nActive), B=propMean[activeOD],
            E=A_pivot[activeLink,activeOD,drop=FALSE], F=Y[tme,activeLink],
            G=diag(nActive), H=rep(0,nActive),
            sdB=propSD[activeOD],
            iter=Xdraws, outputlength=m, burninlength=Xburnin,
            x0=x0Active, type='rda')
        
        if (verbose) {
            cat(sprintf("Accepted ratio = %g\n", X_prop_active$acceptedratio))
        }
        # X_prop <- X_prop$X
        # X_prop <- X_prop[seq(1,Xdraws,floor(Xdraws/m)),]
        
        if (nActive == k) {
            X_prop <- X_prop_active$X
        } else {
            X_prop <- matrix(0, m, k)
            X_prop[,activeOD] <- X_prop_active$X
        }
        
        # Calculate weights
        w <- rep(0,m)

        # Drawing lambda from prior; drawing X from truncated normal on
        # feasible region
        meanMat <- exp(lambda_prop)
        # meanMat <- exp(lambda_prop + phi_prop^2/2)
        # sdMat <- ( meanMat * sqrt(exp(phi_prop^2) - 1) )
        sdMat <- sqrt( meanMat^tau * (exp(phi_prop^2) - 1) )
        truncMat <- matrix(
            pnorm(0, meanMat, sdMat, lower.tail=FALSE, log.p=TRUE),
            m, k)
        
        w <- w + rowSums(matrix(
            dnorm(X_prop, meanMat, sdMat, log=TRUE) - truncMat,
            m, k))
        w <- w - colSums(matrix(
            dnorm( t(X_prop), propMean, propSD, log=TRUE),
            k, m))

        w <- exp(w - max(w))
        # Handle numerical errors
        w[is.na(w)] <- 0
        
        # Normalize importance weights
        w <- w / sum(w)

        # Check effective number of particles
        nEff <- 1/sum(w*w)

        # If verbose, print importance weight diagnostics
        if (verbose) {
            cat(sprintf('Effective number of particles: %g\n', nEff))
            if (nEff < nThresh)
                cat('Redrawing lambda from prior\n')
        }

        # If too few particles, redraw lambda from prior
        if (nEff < nThresh && tme > 1) {
            meanPrior <- colSums( rho^(tme:1 - 1) * prior$mu[1:tme,] )
            sdPrior <- 1/(1-rho^2) * prior$sigma[tme,]

            lambda_prop <- matrix( rnorm(k*m,
                colSums( rho^(tme:1 - 1) * prior$mu[1:tme,] ),
                1/(1-rho^2) * prior$sigma[tme,]),
                m, k, byrow=TRUE )

            # Draw X's from truncated normal on feasible region via xsample, 
            # then thin
            propMean <- exp( meanPrior )
            propSD <- sqrt( (exp(prior$phi$scale[tme]/prior$phi$df)-1) *
                propMean^tau )
                
            # Start with IPF beginning from propMean
            x0 <- ipf(Y[tme,], A_pivot, propMean, tol=1e-6, maxit=Xdraws*1e2,
                verbose=FALSE)
            
            # Now, lsei for refinement
            x0Active <- lsei(
                A=diag(nActive), B=x0[activeOD],
                E=A_pivot[activeLink,activeOD], F=Y[tme,activeLink],
                G=diag(nActive), H=rep(0,nActive),
                tol=sqrt(.Machine$double.eps),
                type=1)$X
            
#            # Start with mirror algorithm
#            # Wrong distribution (uniform, not normal) and slow, but
#            # gets out of corner that lsei found
#            x0Active <- xsample(
#                E=A_pivot[activeLink,activeOD], F=Y[tme,activeLink],
#                G=diag(nActive), H=rep(0,nActive),
#                iter=10,
#                x0=x0Active, type='mirror')$X
#            x0Active <- drop(tail(x0Active, 1))
#            
#            if (verbose) {
#                cat(sprintf("Mirror complete; starting RDA\n"))
#            }
            
            # Switch to RDA algorithm
            # Correct distribution (truncated normal) and faster
            X_prop_active <- xsample(
                A=diag(nActive), B=propMean[activeOD],
                E=A_pivot[activeLink,activeOD], F=Y[tme,activeLink],
                G=diag(nActive), H=rep(0,nActive),
                sdB=propSD[activeOD],
                iter=Xdraws, outputlength=m, burninlength=Xburnin,
                x0=x0Active, type='rda')
            
            if (verbose) {
                cat(sprintf("Accepted ratio = %g\n", X_prop_active$acceptedratio))
            }
            # X_prop <- X_prop$X
            # X_prop <- X_prop[seq(1,Xdraws,floor(Xdraws/m)),]
            
            if (nActive == k) {
                X_prop <- X_prop_active$X
            } else {
                X_prop <- matrix(0, m, k)
                X_prop[,activeOD] <- X_prop_active$X
            }
            
            # Drawing lambda from prior; drawing X from truncated normal on
            # feasible region
            meanMat <- exp(lambda_prop)
            # meanMat <- exp(lambda_prop + phi_prop^2/2)
            # sdMat <- ( meanMat * sqrt(exp(phi_prop^2) - 1) )
            sdMat <- sqrt( meanMat^tau * (exp(phi_prop^2) - 1) )
            truncMat <- matrix(
                pnorm(0, meanMat, sdMat, lower.tail=FALSE, log.p=TRUE),
                m, k)
            
            w <- rep(0, m)
            w <- w + rowSums(matrix(
                dnorm(X_prop, meanMat, sdMat, log=TRUE) - truncMat,
                m, k))
            w <- w - colSums(matrix(
                dnorm( t(X_prop), propMean, propSD, log=TRUE),
                k, m))
            
            w <- exp(w - max(w))
            # Handle numerical errors
            w[is.na(w)] <- 0
            
            # Normalize importance weights
            w <- w / sum(w)

            # Check effective number of particles
            nEff <- 1/sum(w*w)

            if (verbose) {
                cat(sprintf('Effective number of particles after redraw: %g\n',
                    nEff))
            }
        }

        # Resample particles according to weights
        ind <- sample(m, m, replace=TRUE, prob=w)
        
        lambda <- lambda_prop[ind,]
        phi <- phi_prop[ind]
        X <- X_prop[ind,]

        # Move step
        move.out <- move_step(y=Y[tme,], X=X, tme,
            lambda=lambda, phi=phi, lambdatm1=lambdatm1, phitm1=phitm1,
            prior=prior,
            A=A_pivot, A1_inv=A1_inv, A2=A2,
            rho=rho, tau=tau,
            m=m, l=l, k=k,
            ndraws=Movedraws, verbose=verbose)

        # Get revised results
        lambda <- move.out$lambda
        X <- move.out$X
        phi <- move.out$phi

        # Store results from this iteration
        nEffVec[tme] <- nEff
        
        if (aggregate) {
            xList[[tme]] <- agg(X)
            lambdaList[[tme]] <- agg(lambda)
            phiList[[tme]] <- agg(phi)
        } else {
            xList[[tme]] <- X
            lambdaList[[tme]] <- lambda
            phiList[[tme]] <- phi
        }

        # If verbose, print parameter summaries
        if (verbose) {
            endTime <- proc.time()
            print(mean(phi))
            print( rbind(exp(colMeans(lambda)), colMeans(X),
                apply(X,2,median)) )
            cat(sprintf("Runtime for iteration %d:\n", tme))
            print(endTime-startTime)
        }
    }
    
    # Un-reverse results, if needed
    if (backward) {
        Y <- Y[n:1,]
        # prior$mu <- prior$mu[n:1,]
        # prior$sigma <- prior$sigma[n:1,]

        xList <- xList[n:1]
        lambdaList <- lambdaList[n:1]
        phiList <- phiList[n:1]
    }

    # Return list containing results
    unpivot <- numeric(k)
    unpivot[pivot] <- seq(k)

    xList <- lapply(xList, function(mat) mat[,unpivot])
    lambdaList <- lapply(lambdaList, function(mat) mat[,unpivot])

    retval <- list( xList=xList, lambdaList=lambdaList, phiList=phiList,
        y=y, rho=rho, prior=prior, n=n, l=l, k=k,
        A=A, A_qr=A_qr, A1=A1, A1_inv=A1_inv, A2=A2,
        nEff=nEffVec,
        tStart=tStart, backward=backward, aggregate=aggregate )
    
    return(retval)
    # return(list( Xhat=Xhat[,unpivot], Xmed=Xmed[,unpivot],
    #     lambdahat=lambdahat[,unpivot], phihat=phihat ))
}

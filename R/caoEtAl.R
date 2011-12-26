# Functions for EM estimation of network tomography
# Based on methods of Cao et al. (JASA, 2000)


# Check for deterministically-known OD flows at single time
getActive <- function(y, A) {
    odRanges <- xranges(E=A, F=y, ispos=TRUE)
    activeOD <- (odRanges[,2]-odRanges[,1]>0)
    return( activeOD )
}

lambda_init <- function(Y, A, full=FALSE) {
    J <- 2*sqrt(ncol(A))
    
    # Append final row if needed
    if (!full) {
        A <- rbind(A, 2-colSums(A))
    }

    # Aggregate data
    link.means <- colMeans(Y)
    if (!full) {
        link.means <- c(link.means,
                        sum(link.means[1:J/2])-sum(link.means[J/2:ncol(Y)]))
    }

    lambda0 <- rep(0,ncol(A))

    for (i in 1:(J/2)) {
        # Use source data
        lambda0[seq((i-1)*J/2+1, i*J/2)] <- lambda0[seq((i-1)*J/2+1, i*J/2)] + 
            link.means[i]/(J)

        # Use destination data
        lambda0[seq(i,ncol(A),J/2)] <- ( lambda0[seq(i,ncol(A),
            J/2)] + link.means[i+J/2]/(J) )
    }

    return(lambda0)
}

phi_init <- function(Y, A, lambda0, c) {
    exp(mean(log(diag(cov(Y))/ (A %*% lambda0^c) )))
}

Q_iid <- function(logtheta, c, M, rdiag, epsilon) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])
    phi <- exp(logtheta[length(logtheta)])

    # Calculate first part of Q :
    #  -1/2 sum_t (m_t - lambda)' \Sigma^{-1} (m_t - lambda)
    Q <- -1/2*sum(apply( M, 1, function(m) sum( (m-lambda)^2 /
        (phi*lambda^c+epsilon) ) ))

    # Calculate remaining part:
    ## -T/2 (I*log(\phi) + \sum_i c*log(\lambda_i) + 
    # \sum_i r_ii/ \phi / lambda_i^c)
    # With nugget, revised to:
    # -T/2 ( \sum_i log(\phi * \lambda_i^c + \varepsilon) +
    # \sum_i r_ii/(\phi * \lambda_i^c + \varepsilon) )
    Q <- Q - nrow(M)/2*( sum(log( phi*lambda^c + epsilon )) +
        sum(rdiag/(phi * lambda^c + epsilon) ) )

    return(Q);
}

m_estep <- function(yt, lambda, phi, A, c, epsilon) {
    sigma <- diag_mat(phi*lambda^c + epsilon)
    
    m <- lambda
    m <- m  +  sigma %*% t(A) %*%
        # solve(A %*% diag(phi*lambda^c) %*% t(A)) %*%
        # chol2inv(chol(A %*% diag(phi*lambda^c + epsilon) %*% t(A))) %*%
        solve(A %*% sigma %*% t(A), yt - A %*% lambda)
    
    return(m)
}

R_estep <- function(lambda, phi, A, c, epsilon) {
    sigma <- diag_mat(phi*lambda^c + epsilon)
    AdotSigma <- A %*% sigma
    
    R <- diag(phi*lambda^c + epsilon)
    R <- R - t(AdotSigma) %*%
        # solve(A %*% diag(phi*lambda^c) %*% t(A)) %*% A %*%
        # chol2inv(chol(A %*% diag(phi*lambda^c + epsilon) %*% t(A))) %*% )
        solve( AdotSigma %*% t(A), AdotSigma )
    
    return(R)
}

grad_iid <- function(logtheta, c, M, rdiag, epsilon) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])
    phi <- exp(logtheta[length(logtheta)])

    grad <- rep(NA, length(logtheta))

    # Calculate gradient for phi
    grad[length(grad)] <- ( -nrow(M)/2*(ncol(M)/phi -
        1/phi/phi*sum(rdiag/lambda^c)) + 1/2/phi/phi*
        sum(apply(M, 1, function(m) sum((m-lambda)^2/lambda^c))) )

    # Calculate gradient for each lambda
    grad[-length(grad)] <- (-nrow(M)/2*c*(1/lambda - rdiag/phi/lambda^(c+1))
        + (colSums(M)-nrow(M)*lambda)/phi/lambda^c
        + c/2/phi*sapply(1:ncol(M),
        function(j) sum((M[,j]-lambda[j])^2/lambda[j]^(c+1))))

    # Jacobian adjustment
    grad <- grad*exp(logtheta)

    # Return gradient
    return(grad)
}

locally_iid_EM <- function(Y, c=2, A=NULL, lambda0=NULL, phi0=NULL,
    maxiter = 1e3, tol=1e-6, epsilon=0.01, full=FALSE) {
    
    # Check for inactive (deterministically-known) OD flows
    activeLink <- which(apply(Y, 2, function(col) max(col) > 0))
    activeOD <- apply(Y, 1, getActive, A=A)
    activeOD <- which(apply(activeOD, 1, any))
    
    # Determine number of latent variables
    p <- 1
    
    if (is.null(A)) {
        p <- (ncol(Y)/2)^2
        A <- make_routemat(p)
        A <- A[-nrow(A),]
    } else {
        p <- ncol(A)
    }

    # Setup initial parameters
    if (is.null(lambda0)) {
        lambda0 <- lambda_init(Y,A,full=full)
        lambda0[!is.finite(lambda0)] <- 1
        lambda0[lambda0<=0] <- 1
        lambda <- lambda0
    } else {
        lambda <- lambda0
    }
    
    # Subset A and lambda
    A <- A[activeLink,activeOD]
    Y <- Y[,activeLink]
    lambda0 <- lambda0[activeOD]
    lambda <- lambda[activeOD]
    
    if (is.null(phi0)) {
        phi0 <- phi_init(Y, A, lambda0, c)
        phi <- phi0
    } else {
        phi <- phi0
    }
    
    # Calculate bounds
    lower <- rep(log(.Machine$double.eps)/c+log(max(Y)), p)
    upper <- rep(log(max(Y)), p)

    # Run first iteration
    # E step - calculation of m_t = E(x_t | y_t, \theta) and 
    #  R = Var(x_t | y_t, \theta)
    tryCatch(M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi, A=A, c=c,
        epsilon=epsilon)),
        error = function(e) {print(phi); print(lambda); print("M"); stop(e);},
        finally = NULL )
    R <- tryCatch( R_estep(lambda, phi=phi, A, c, epsilon),
        error = function(e) {print(phi); print(lambda); print("M"); stop(e);},
        finally = NULL )
    # M step - Optimize Q wrt theta
    mstep <- optim(log(c(lambda0,phi0)),
        function(...) -Q_iid(...),# gr=function(...) -grad_iid(...),
        c=c, M=M, rdiag=diag(R), epsilon=epsilon,
        lower=lower, upper=upper,
        method="L-BFGS-B" )
    theta <- exp(mstep$par)
    lambda <- theta[-length(theta)]
    phi <- theta[length(theta)]

    # Calculate initial ll
    ll <- -mstep$value
    
    # Run EM iterations
    for (iter in 1:maxiter) {
        # E step - calculation of m_t = E(x_t | y_t, \theta) and 
        #  R = Var(x_t | y_t, \theta)
        tryCatch(M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi, A=A, c=c,
        epsilon=epsilon)),
        error = function(e) {print(phi); print(lambda); print("M"); stop(e);},
        finally = NULL )
        R <- tryCatch( R_estep(lambda, phi=phi, A, c, epsilon),
            error = function(e) {print(phi); print(lambda); print("M"); stop(e);},
            finally = NULL )
        
        # M step - Optimize Q wrt theta
        mstep <- optim(log(c(lambda0,phi0)),
            function(...) -Q_iid(...), gr=function(...) -grad_iid(...),
            c=c, M=M, rdiag=diag(R), epsilon=epsilon,
            lower=lower, upper=upper,
            method="L-BFGS-B" )
        theta <- exp(mstep$par)
        lambda <- theta[-length(theta)]
        phi <- theta[length(theta)]
        ll.new <- -mstep$value

        # cat(ll, ll.new, abs(ll.new-ll)/abs(ll.new+ll)*2, "\n")

        # Check for convergence
        if (abs(ll.new-ll)/abs(ll.new+ll)*2 < tol) {
            ll <- ll.new
            break
        } else {
            ll <- ll.new
        }
    }
    tmp <- numeric(p)
    tmp[ activeOD ] <- lambda
    lambda <- tmp
    return(list(lambda=lambda, phi=phi, iter=iter))
}

Q_smoothed <- function(logtheta, c, M, rdiag, eta0, sigma0, V,
    eps.lambda, eps.phi) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])+eps.lambda
    phi <- exp(logtheta[length(logtheta)])+eps.phi
    logtheta <- c(log(lambda),log(phi))
    
    # Calculate first part of Q (from locally iid model)
    Q <- Q_iid(logtheta, c, M, rdiag)

    # Add log-normal prior component
    Q <- (Q - 1/2*t(logtheta-eta0) %*% chol2inv(chol(V+sigma0))
        %*% (logtheta-eta0))

    return(Q);
}

grad_smoothed <- function(logtheta, c, M, rdiag, eta0, sigma0, V,
    eps.lambda, eps.phi) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])+eps.lambda
    phi <- exp(logtheta[length(logtheta)])+eps.phi
    logtheta <- c(log(lambda),log(phi))

    # Calculate first part of grad (from iid model)
    grad <- grad_iid(logtheta, c, M, rdiag)

    # Add gradient from normal prior
    grad <- grad - chol2inv(chol(V+sigma0)) %*% (logtheta - eta0)

    # Return gradient
    return(grad)
}


smoothed_EM <- function(data=Y, eta0, sigma0, V, c=2, A=NULL,
    maxiter = 1e3, tol=1e-6) {
    eps.lambda <- 0
    eps.phi <- 0
    # Determine number of latent variables
    p <- (ncol(Y)/2)^2
    
    if (is.null(A)) {
        A <- make_routemat(p)
        A <- A[-nrow(A),]
    }

    # Setup initial parameters
    theta0 <- exp(eta0)

    lambda0 <- theta0[-length(eta0)] + 1e-2
    # lambda0 <- (lambda0 + lambda_init(Y))/2
    # lambda0 <- lambda_init(Y)
    lambda <- lambda0

    phi0 <- theta0[length(eta0)] + 1e-2
    phi <- phi0

    # Remove final column from Y to avoid colinearity
    Y <- as.matrix(Y[,-ncol(Y)])

    # Run first iteration
    # E step - calculation of m_t = E(x_t | y_t, \theta) and 
    #  R = Var(x_t | y_t, \theta)
    M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi,
        A=A, c=c))
    R <- R_estep(lambda, phi=phi, A, c)
    
    # M step - Optimize Q wrt theta
    mstep <- optim(c(log(lambda0),log(phi0)),
        function(...) -Q_smoothed(...),
        gr=function(...) -grad_smoothed(...),
        c=c, M=M, rdiag=diag(R), eta0=eta0, sigma0=sigma0, V=V,
        eps.lambda=eps.lambda, eps.phi=eps.phi,
        method="BFGS")
    theta <- exp(mstep$par)
    lambda <- theta[-length(theta)]+eps.lambda
    phi <- theta[length(theta)]+eps.phi

    # Calculate initial ll
    ll <- -mstep$value
    
    # Run EM iterations
    for (iter in 1:maxiter) {
        # E step - calculation of m_t = E(x_t | y_t, \theta) and 
        #  R = Var(x_t | y_t, \theta)
        M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi,
            A=A, c=c))
        R <- R_estep(lambda, phi=phi, A, c)
        
        # M step - Optimize Q wrt theta
        mstep <- optim(c(log(lambda0),log(phi0)),
            function(...) -Q_smoothed(...),
            gr=function(...) -grad_smoothed(...),
            c=c, M=M, rdiag=diag(R), eta0=eta0, sigma0=sigma0, V=V,
            eps.lambda=eps.lambda, eps.phi=eps.phi,
            method="BFGS")
        theta <- exp(mstep$par)
        lambda <- theta[-length(theta)]+eps.lambda
        phi <- theta[length(theta)]+eps.phi
        ll.new <- -mstep$value

        # cat(ll, ll.new, abs(ll.new-ll)/abs(ll.new+ll)*2, "\n")

        # Check for convergence
        if (abs(ll.new-ll)/abs(ll.new+ll)*2 < tol) {
            ll <- ll.new
            break
        } else {
            ll <- ll.new
        }
    }
    mstep <- optim(log(theta),
        function(...) -Q_smoothed(...),
        gr=function(...) -grad_smoothed(...),
        c=c, M=M, rdiag=diag(R), eta0=eta0, sigma0=sigma0, V=V,
        eps.lambda=eps.lambda, eps.phi=eps.phi,
        method="BFGS", hessian=TRUE)
    
    return(list(lambda=lambda, phi=phi, iter=iter, etat=c(log(lambda),log(phi)),
        sigmat=chol2inv(chol(mstep$hessian))))
}

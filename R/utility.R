#' Convert time string to decimal hour
#'
#' @param x input character vector of times
#' @param fmt input character format for times
#' @return numeric vector of decimal times in hours
#' @keywords character
#' @export
#' @examples
#' strphour("31/08/87 12:53:29")
#'
strphour <- function(x, fmt="(%m/%d/%y %H:%M:%S)") {
    ptime <- strptime(as.character(x), format=fmt)
    ptime$hour + ptime$min/60 + ptime$sec/3600
}

# Build routing matrices for 1router & 2router data
make_routemat <- function(p) {
    # Allocate matrix
    J <- sqrt(p)*2
    routemat <- matrix(0, J, p)
    
    # Setup nonzero entries
    for (i in 1:(nrow(routemat)/2)) {
        routemat[i,seq((i-1)*J/2+1, i*J/2)] <- 1
        routemat[i+nrow(routemat)/2, seq(i,ncol(routemat), J/2)] <- 1
    }

    return(routemat)
}

# Functions to handle diagonal matrices much faster (array is very slow)
# Make diagonal matrix from vector
diag_mat <- function(x) {
    n <- length(x)
    y <- matrix(0,n,n)
    y[1L + 0L:(n-1L) * (n+1L)] <- x
    return(y)
}

# Make vector of 1-dimensional diagonal indices
diag_ind <- function(n) {
    1L + 0L:(n-1L)*(n+1L)
}

# Thin vector of indices for MCMC
thin <- function(m, interval=10) {
    seq(1,m,interval)
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


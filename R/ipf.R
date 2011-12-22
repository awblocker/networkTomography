# Load shared libraries
dyn.load("src/ipf.so")

ipf <- function(y, A, x0,
    tol=.Machine$double.eps, maxit=1e3, verbose=FALSE) {
    # Get active rows
    activeRows <- which(y > 0)
    
    # Zero inactive columns
    if ( any(y==0) ) {
        activeCols <- !pmin(1, colSums(A[y==0,,drop=FALSE]))
    } else {
        activeCols <- rep(TRUE, ncol(A))
    }
    x0[!activeCols] <- 0
    x0[activeCols] <- pmax(1, x0[activeCols])
    
    # Run IPF
    ans <- .Call("ipf", y[activeRows], A[activeRows, activeCols, drop=FALSE],
            dim(A[activeRows, activeCols, drop=FALSE]), x0[activeCols],
            as.numeric(tol), as.integer(maxit), as.logical(verbose))
    
    x0[activeCols] <- ans$x
    
    return(x0)
}

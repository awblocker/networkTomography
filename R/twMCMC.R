
# Function to run MCMC algorithm of Tebaldi & West 1998
network_mcmc <- function(Y, A, prior, ndraws=1.2e5, burnin=2e4, verbose=0,
                         obj_param=3) {
    # Format verification
    Y <- as.numeric(Y)

    # Initialize X and lambda
    n <- ncol(A)
    k <- nrow(A)
    X <- matrix(0, ndraws+burnin, n)

    # Decompose A matrix into full-rank and singular parts; retain pivot info
    A_qr <- qr(A)
    pivot <- A_qr$pivot
    A_pivot <- A[,pivot]
    A1 <- A_pivot[,seq(A_qr$rank)]
    A1_inv <- solve(A1)
    A2 <- A_pivot[,seq(A_qr$rank+1,n)]

    # Pivot prior
    prior$a <- prior$a[pivot]
    prior$b <- prior$b[pivot]

    # Initialize X via simplex method (gives integer solutions with linear
    # constraints for transportation problems)
    # obj <- rpois(n, obj_param)*(rbinom(n,1,1/2)*2 - 1)
    obj <- rnorm(n, 0, obj_param)
    const_mat <- rbind(A_pivot, diag(n))
    const_rhs <- c(Y, rep(1,n))
    const_dir <- c(rep('==',k), rep('>=',n))
    int_vec <- seq(n)

    #init <- lp(objective.in=obj, const.mat=const_mat, const.dir=const_dir,
    #    const.rhs=const_rhs, int.vec=int_vec)
    init <- Rglpk_solve_LP(obj, const_mat, const_dir, const_rhs)
    X[1,] <- init$solution
    X1 <- X[1,seq(k)]
    X2 <- X[1,seq(k+1,n)]

    # If verbose, print initial solution
    if (verbose) print(X[1,])

    # Initialize lambda with draw from conditional posterior
    lambda <- matrix(0, ndraws+burnin, n)
    lambda[1,] <- rgamma(n, prior$a+X[1,], prior$b+1)

    # Setup adjList
    adjList <- lapply(1:(n-k), function(j) A1_inv %*% A2[,j])
    posList <- lapply(adjList, function(x) which(x>0))
    negList <- lapply(adjList, function(x) which(x<0))

    # Loop for MCMC iterations
    accepts <- rep(0, n-k)
    for (iter in seq(2,ndraws+burnin)) {
        # Draw X2_j, lambda_j+k
        for (j in seq(n-k)) {
            # MH step on X2_j

            # Propose uniformly along feasible region
            adjVec <- adjList[[j]]
            posVec <- posList[[j]]
            negVec <- negList[[j]]
            remainderVec <- X1 + adjVec * X2[j]
            limitVec <- remainderVec/adjVec
            maxX2 <- max(min(limitVec[posVec]), 0)
            minX2 <- max(max(limitVec[negVec]), 0)
            proposal <- sample(floor(maxX2)-ceiling(minX2)+1, 1)
            proposal <- proposal - 1 + ceiling(minX2)

            # Propose lambda_j+k | X2_j
            lambdaProp <- rgamma(1, proposal + prior$a[j+k] + 1,
                                 prior$b[j+k] + 1)

            # If feasible, MH step
            llr <- (log(lambdaProp)*proposal - lgamma(proposal+1) -
                    log(lambda[iter-1,j+k])*X2[j] + lgamma(X2[j]+1) )
            # Check feasibility (issues with integer constraints)
            llr <- llr + log(min(remainderVec - adjVec*proposal)>0)
            lir <- 0

            # Accept/reject part
            if (!is.na(llr-lir)) {
                if (log(runif(1)) < llr - lir) {
                    X2[j] <- proposal
                    accepts[j] <- accepts[j] + 1
                    X1 <- remainderVec - adjVec * X2[j]
                    lambda[iter,j+k] <- lambdaProp
                }
            }
        }

        X[iter,] <- c(X1,X2)

        # Draw lambda_1, ..., lambda_k
        lambda[iter,1:k] <- rgamma(k, X[iter,1:k] + prior$a[1:k] + 1,
                                   prior$b[1:k] + 1)

        # If verbose, print updates every 1e3 iterations
        if (verbose>0 & iter %% 1e3 == 0) {
            cat(sprintf('Iter %d\n', iter))
            if (verbose > 1) {
                print(lambda[iter,])
                print(X[iter,])
            }
        }
    }

    inv <- numeric(n)
    inv[pivot] <- seq(n)

    lambda <- lambda[,inv]
    X <- X[,inv]
    accepts <- accepts[inv]

    out <- tail(data.frame(X, lambda),ndraws)
    colnames(out) <- c(paste('X', seq(n), sep=''),
                       paste('lambda', seq(n), sep='') )

    return(list(lambda_out=mcmc(out[,seq(n+1,2*n)]), X_out=mcmc(out[,seq(n)]),
                accepts=accepts))
}


vardi.algorithm <-
function(A, y, lambda, B, S) {
  epsilon <- 1e-3
  lambda.new <- vardi.iteration(A,y,lambda,B,S)
  while (!identical(abs(lambda.new - lambda) > epsilon,
                    rep(FALSE,length(lambda)))) {
    lambda <- lambda.new
    lambda.new <- vardi.iteration(A,y,lambda,B,S)
  }

  return(lambda.new)
}


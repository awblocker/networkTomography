vardi.compute.BS <- function(A,Y) {
  r <- dim(A)[1]
  K <- dim(Y)[2]
  first <- 1
  for (i in 1:r) {
    for (j in i:r) {
      bi <- A[i,]*A[j,]
      si <- 1/K * sum(Y[i,]*Y[j,]) - mean(Y[i,])*mean(Y[j,])
      # Delete the impossible moment equations under the Poisson model.
      if (si < 0 && sum(bi) > 0 || si > 0 && sum(bi) < 0 ||
          si != 0 && sum(bi) == 0 || si == 0 && sum(bi) != 0) {
        next
      }
      if (first) {
        B <- bi
        S <- si
        first <- 0
        next
      }
      B <- rbind(B,bi)
      S <- c(S,si)
    }
  }
  dimnames(B) <- NULL

  return(list(B,S))
}

vardi.iteration <- function(A, y, lambda, B, S) {
  adot <- apply(A,2,sum)
  bdot <- apply(B,2,sum)
  Alambda <- as.vector(crossprod(t(A),lambda))    # A %*% lambda
  Blambda <- as.vector(crossprod(t(B),lambda))    # B %*% lambda
  Ay <- A*y
  BS <- B*S
  AyoAlambda <- Ay/Alambda
  BSoBlambda <- BS/Blambda
  lambda.hat.Ayl <- lambda/apply(A,2,sum)*apply(AyoAlambda,2,sum)
  lambda.hat.BSl <- lambda/apply(B,2,sum)*apply(BSoBlambda,2,sum)

  lambda.new <- adot/(adot+bdot)*lambda.hat.Ayl +
    bdot/(adot+bdot)*lambda.hat.BSl

  return(lambda.new)
}

vardi.algorithm <- function(A, y, lambda, B, S) {
  epsilon <- 1e-3
  lambda.new <- vardi.iteration(A,y,lambda,B,S)
  while (!identical(abs(lambda.new - lambda) > epsilon,
                    rep(FALSE,length(lambda)))) {
    lambda <- lambda.new
    lambda.new <- vardi.iteration(A,y,lambda,B,S)
  }

  return(lambda.new)
}


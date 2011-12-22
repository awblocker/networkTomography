vardi.iteration <-
function(A, y, lambda, B, S) {
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


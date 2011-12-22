vardi.compute.BS <-
function(A,Y) {
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


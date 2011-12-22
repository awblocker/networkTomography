\name{compute.BS}
\alias{compute.BS}
\title{
Function that takes input matricies A and Y from the example1.R file and calculates matricies B and S
}
\description{
This is a middleman step in the algorithm that helps clean up the rest of the code.
}
\usage{
compute.BS(A, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

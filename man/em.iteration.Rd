\name{em.iteration}
\alias{em.iteration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that approximates lambda using method of moments estimation as described in Vardi 1996. This function must be iterated (which is acomplished with em.algorithm.R) in order to have the ouputs converge to the actual estimation, lambda hat. For examples, see Vardi.Algorithm.Rd
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
em.iteration(A, y, lambda, B, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
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
function(A, y, lambda, B, S) {
  adot <- apply(A,2,sum)
  bdot <- apply(B,2,sum)
  Alambda <- as.vector(crossprod(t(A),lambda))    # A \%*\% lambda
  Blambda <- as.vector(crossprod(t(B),lambda))    # B \%*\% lambda
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{Y}
\alias{Y}
\docType{data}
\title{
Traffic intensity input across the network's direct links
}
\description{
This input matrix, who's columns correspond to the time divisions and who's rows represent direct links, is one of the two acutal inputs (the other being matrix A) that are needed to utilize the algorithm rather than test it. Y = AX.
}
\usage{data(Y)}
\format{
  The format is:
 num [1:7, 1:100] 1 7 30 6 47 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Y)
## maybe str(Y) ; plot(Y) ...
}
\keyword{datasets}

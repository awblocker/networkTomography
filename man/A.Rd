\name{A}
\alias{A}
\docType{data}
\title{
Network Matrix
}
\description{
r x c matrix where r is the number of direct links and c is the number of Object Destination pairs. Given n nodes in a network, c is generally equal to n(n-1). r is dependent on the particular matrix formation. A number greater than zero indicates that for a message to travel along the corresponding column's path it could (or must if the value =1) travel through the corresponding row's direct link.
}
\usage{data(A)}
\format{
  The format is:
 num [1:7, 1:12] 1 0 0 0 0 0 0 0 1 0 ...
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
data(A)
## maybe str(A) ; plot(A) ...
}
\keyword{datasets}

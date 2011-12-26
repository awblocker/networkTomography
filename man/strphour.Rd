\name{strphour}
\alias{strphour}
\title{Convert time string to decimal hour}
\usage{
  strphour(x, fmt = "(\%m/\%d/\%y \%H:\%M:\%S)")
}
\arguments{
  \item{x}{input character vector of times}

  \item{fmt}{input character format for times}
}
\value{
  numeric vector of decimal times in hours
}
\description{
  Convert time string to decimal hour
}
\examples{
strphour("31/08/87 12:53:29")
}
\keyword{character}


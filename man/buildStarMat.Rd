\name{buildStarMat}
\alias{buildStarMat}
\title{Build routing matrix for star network topology}
\usage{
buildStarMat(n)
}
\arguments{
  \item{n}{integer number of nodes in the network}
}
\value{
matrix of dimension 2n x n^2 that transforms OD flows to
link loads
}
\description{
Build routing matrix for star network topology
}
\examples{
buildStarMat(3)
}
\keyword{array}


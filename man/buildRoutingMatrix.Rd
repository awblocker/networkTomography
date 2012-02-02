\name{buildRoutingMatrix}
\alias{buildRoutingMatrix}
\title{Build routing matrix from table of link relationships}
\usage{
  buildRoutingMatrix(nodes, src, dest, weights = NULL,
    agg = FALSE, sep = "_", aggChar = "*", verbose = 0)
}
\arguments{
  \item{nodes}{vector (lenght n) of node identifiers}

  \item{src}{vector (length m) of sources, one per link,
  matched with dest}

  \item{dest}{vector (length n) of destination identifiers,
  one per link, matched with src}

  \item{weights}{numeric vector (length m) of weights for
  each link; used in shortest-path routing calculations
  (roughly OSPF)}

  \item{agg}{logical for whether to include aggregate
  source and destination flows for each node}

  \item{sep}{character separator between node id's for link
  and OD names}

  \item{aggChar}{character to indicate aggregate flows;
  should be distinct from sep}

  \item{verbose}{integer level of verbosity; 0 is silent,
  >=1 are increasing levels of reporting}
}
\value{
  List consisting of routing matrix \code{A} (dense) of
  dimensions m x n and iGraph object for network
  \code{topo}
}
\description{
  Constructs routing matrix from link relationships.
  Determines routes using (weighted) shortest-path
  calculation (mirroring OSPF). Currently handles tied
  paths arbitrarily; will incorporate fractions for tie
  resolution in next version. Can optionally include
  aggregate source and destination flows for each node;
  this can make a major difference for some topologies.
  Tomogravity methods typically make use of such
  information, which most routers collect. Note that
  resulting routing matrix need not be of full row rank.
}
\keyword{array}


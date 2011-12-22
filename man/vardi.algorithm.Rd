\name{em.algorithm}
\alias{em.algorithm}
\title{
function that iterates em.iteration.R to determine lambda hat
}
\description{
This is the last step in the algorithm, which uses the asymptotic nature of em.iteration.R to repeat the iteration until each subsequent value derived from each iteration is less that .001 different from the previous. It then outputs the new lambda values.
}
\usage{
em.algorithm(A, y, lambda, B, S)
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
# Example 1.1

# The routing matrix A is given by:
A <- rbind(c(1,0,0,0,0,0,0,0,0,0,0,0),
           c(0,1,1,0,0,1,0,0,0,0,0,0),
           c(0,0,0,1,0,1,1,0,0,1,0,0),
           c(0,0,0,0,1,0,0,0,0,0,0,0),
           c(0,0,0,0,0,0,1,1,0,1,1,0),
           c(0,0,1,0,0,1,0,0,1,0,0,0),
           c(0,0,0,0,0,0,0,0,0,1,1,1))

# Repeat this whole process 50 times to produce 50 such estimates:
X <- matrix(0,12,100)
Y <- matrix(0,7,100)
Lambda <- matrix(0,12,50)
for (i in 1:50) {
  # "Daily" transmission rate between all 12 pairs:
  lambda <- 1:12

  # Daily data on the seven links for 100 days:
  for (k in 1:100) {
    X[,k] <- rpois(12, lambda)
    Y[,k] <- A %*% X[,k]
  }
  # Make all entries in Y positive:
  Y <- Y + .Machine$double.eps^0.5

  # Calculate the matrix B and covariances S:
  BS <- vardi.compute.BS(A,Y)
  B <- BS[[1]]
  S <- BS[[2]]

  # Apply the Vardi algorithm:
  Lambda[,i] <- vardi.algorithm(A,rowMeans(Y),lambda,B,S)
}

# Example with 1Router Data

library("lattice")
source("~/Dropbox/networkTomographyPkgShare/networkTomography/R/Vardi.algorithm.R")
source("~/Dropbox/networkTomographyPkgShare/networkTomography/R/vardi.compute.BS.R")
source("~/Dropbox/networkTomographyPkgShare/networkTomography/R/vardi.iteration.R")

# Read the complete data from the subnetwork around Router1 in a local network
# at Lucent at regular 5-minute intervals.
router1.data <- read.csv(file="~/Dropbox/networkTomographyPkgShare/networkTomography/data/1router_clean.txt", header=TRUE, row.names=1)
router1.data <- router1.data[,1]
T <- length(router1.data)/25
times <- seq(0,(T-1)*5/60,5/60)
link.loads <- matrix(0,8,T)
od.flows <- matrix(0,16,T)
for (t in 1:T) {
  shift <- 25*(t-1)
  od.flows[,t] <- router1.data[(shift+1):(shift+16)]
  link.loads[,t] <- router1.data[(shift+17):(shift+24)]
}
Y <- link.loads
# Make all entries in Y positive:

# Routing matrix for the single-router network around Router1:
A <- matrix(0,8,16)
for (i in 1:4) {
  A[i,(4*(i-1)+1):(4*i)] <- rep(1,4)
}
for (i in 5:8) {
  A[i,c(i-4,4+i-4,8+i-4,12+i-4)] <- 1
}
Ac <- A
# Discard the last row to make the matrix full rank:
A <- A[1:7,]
J <- dim(A)[1]
I <- dim(A)[2]

# Choose the power c, window size w, and half-width h:
c <- 2
w <- 11
h <- (w-1)/2

# Use the method of moments in the local iid model using a moving window of a
# fixed size w = 2h + 1, where h is the half-width.
ym <- rep(0, 9*T)
Lambda <- matrix(NA,16,T)
for (t in 1:T) {
  start <- t-h
  end <- t+h
  if (start < 1 ) {
    start <- 1
  }
  if (end > T) {
    end <- T
  }
  ym[((t-1)*9+1):((t-1)*9+8)] <- rowMeans(Y[,start:end])
  ym[(t-1)*9+9] <- sum(ym[((t-1)*9+1):((t-1)*9+4)])
  
  # Calculate the matrix B and covariances S:
  BS <- vardi.compute.BS(A,Y[1:7,start:end])
  B <- BS[[1]]
  S <- BS[[2]]

  # Apply the EM algorithm:
  # Choose an initial setting for the parameters lambda
  a0 <- (rep(1,J)  %*% c(apply(Y[1:7,start:end]/T,1,sum))) /
    (rep(1,J) %*% A %*% rep(1,I))
  a0 <- as.numeric(a0)
  lambda <- a0 * rep(1,I)
  Lambda[,t] <- vardi.algorithm(A,rowMeans(Y[1:7,start:end]),lambda,B,S)
}

# Figure 5. Mean Traffic Estimates, \hat{lambda}_t, for all OD Pairs.
timesp <- rep(times[1], 25)
traffic.estimates <- Ac %*% Lambda[1:16,1]
y1 <- c(Lambda[1:16,1], traffic.estimates, sum(traffic.estimates[1:4]))
for (t in 2:T) {
  timesp <- c(timesp, rep(times[t],25))
  traffic.estimates <- Ac %*% Lambda[1:16,t]
  y1 <- c(y1, c(Lambda[1:16,t], traffic.estimates, sum(traffic.estimates[1:4])))
}
for (t in 1:T) {
  timesp <- c(timesp, rep(times[t],9))
}
y1 <- c(y1, ym)
d <- data.frame(x=timesp, y=y1, nme=c(rep(c("fddi->fddi","fddi->switch",
                                  "fddi->local","fddi->corp","switch->fddi",
                                  "switch->switch","switch->local",
                                  "switch->corp","local->fddi","local->switch",
                                  "local->local","local->corp","corp->fddi",
                                  "corp->switch","corp->local","corp->corp",
                                  "origin fddi","origin switch","origin local",
                                  "origin corp","destination fddi",
                                  "destination switch","destination local",
                                  "destination corp","total"),T),
                                  rep(c("origin fddi","origin switch",
                                        "origin local","origin corp",
                                        "destination fddi",
                                        "destination switch",
                                        "destination local","destination corp",
                                        "total"),T)),
                method=c(rep("hat{lambda}_t",25*T),
                  rep("smoothed observation",9*T)))
tp1.traffic.estimates <- xyplot(y~x|nme, data=d, groups=method, type="l",
                                ref=TRUE, auto.key=list(columns=2),
                                aspect=1/1.5, layout=c(5,5,1),
                                between=list(x=c(0,0,0,0.5), y=c(0,0,0,0.5)),
                                xlab="hour of day", xlim=c(0,24),
                                ylab="bytes/sec")
tp1.traffic.estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

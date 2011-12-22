\name{network_mcmc}
\alias{network_mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to run MCMC algorithm of Tebaldi & West 1998
}
\description{
Approaches network tomography problem with Gibbs and Metropolis-Hastings Algorithms and a gamma prior on the poisson liklihood. See details section for more.
}
\usage{
em.iteration(A, y, lambda, B, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Network Matrix as defined by Vardi
}
  \item{y}{
direct link traffic vector
}
  \item{lambda}{
parameter values for each poisson OD pair distribution
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
}
\details{
An iterative process that estimates the OD traffic rates through the following general process: first, initialize X by decomposing the Y = AX equation through QR factorization and then solving for X with the inputted A and Y. Then, use Gibbs sampling techniques, draw Lambda values from the most recent X. Use this new Lambda vector to get a new X vector. Iterate specified number of times (ndraws input)
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
# setting up required packages and functions
library(plyr)
library(Rglpk)
library(slam)
library(coda)
library(foreach)

source("~/Dropbox/networkTomographyPkgShare/networkTomography/R/twMCMC_functions_joint.R")

# Convert time string to decimal hour (needed for loading 1router and 2router
# data)
strphour <- function(x, fmt="(%m/%d/%y %H:%M:%S)") {
    ptime <- strptime(as.character(x), format=fmt)
    ptime$hour + ptime$min/60 + ptime$sec/3600
}

# Build routing matrices for 1router & 2router data
make_routemat <- function(p) {
    # Allocate matrix
    J <- sqrt(p)*2
    routemat <- matrix(0, J, p)
    
    # Setup nonzero entries
    for (i in 1:(nrow(routemat)/2)) {
        routemat[i,seq((i-1)*J/2+1, i*J/2)] <- 1
        routemat[i+nrow(routemat)/2,seq(i,ncol(routemat),J/2)] <- 1
    }

    return(routemat)
}

# Set options
options(cores=2)

# Set parameters
t0 <- 5
ndraws <- 3e4
burnin <- 3e4
nchains <- 2
#nWorkers <- 8
thinInterval <- 100
logOutput <- TRUE

set.seed(t0+20081025)

# Load data
onerouter <- read.table("~/Dropbox/networkTomographyPkgShare/networkTomography/data/1router_clean.txt", header=TRUE, sep=',')
onerouter$time <- sapply(onerouter$time, strphour)
od <- onerouter[onerouter$od==1,]
link <- onerouter[onerouter$od==0,]
link <- link[link$nme!="total",]


# Setup groundtruth OD flow data
groundtruth <- daply(od, .(time), function(df) df$value)

# Setup y and attribute data
y <- daply(link, .(time), function(df) df$value)
n <- nrow(y)
tvec <- daply(link, .(time), function(df) as.numeric(unique(df$time)))

# Rescale and reshape data
y <- y[,-ncol(y)]

k <- ncol(groundtruth)

# Setup routing matrix
routemat <- make_routemat(k)
A <- routemat[-nrow(routemat),]

# set prior
prior <- list(a=rep(0,k),b=rep(0,k))

X_results <- vector("list", nchains)
lambda_results <- vector("list", nchains)
accept_results <- vector("list", nchains)

for (i in seq(nchains)) {
    tStart <- proc.time()
    
    # Run MCMC
    result <- network_mcmc( y[t0,], A, prior, ndraws, burnin, 1 )
   
    #
    X_results[[i]] <- result$X_out
    lambda_results[[i]] <- result$lambda_out
    accept_results[[i]] <- result$accepts
    
    tEnd <- proc.time()
    cat(sprintf("Runtime for chain %d:\n", i), file=stderr())
    cat(tEnd-tStart, "\n", file=stderr())
}

View(X_results)

# visual data representation
# for (i in seq(k)) {
    # hist(sapply(X_results[,i,drop=FALSE],identity), 'fd',
    # xlab=varnames(X_results)[i], ylab='', main='', freq=FALSE)
}
# Example 1.1 for TW

# get TW function
source("~/Dropbox/networkTomographyPkgShare/networkTomography/R/twMCMC_functions_joint.R")

# Set parameters
t0 <- 5
ndraws <- 3e3
burnin <- 2e3
nchains <- 1
#nWorkers <- 8
thinInterval <- 100
logOutput <- TRUE

set.seed(t0+20081025)

k<-100

# set prior
prior <- list(a=rep(0,k),b=rep(1,k))

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
  Y <- Y + .Machine$double.eps^0.5
}

X_results <- vector("list", nchains)
lambda_results <- vector("list", nchains)
accept_results <- vector("list", nchains)

# Run MCMC
result <- network_mcmc( Y[t0,], A, prior, ndraws, burnin, 2 )
   
# get Lambda result
LambdaMeans<-rep(0,12)
for (i in 1:12) {
LambdaMeans[i]<-mean(result$lambda_out[,i])
}
LambdaMeans

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

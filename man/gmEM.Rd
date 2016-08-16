\name{gmEM}
\alias{gmEM}
\title{
Expectation-Maximization Algorithm
}
\description{
Function for deriving a Markov generator matrix estimate by an instance of the expectation-maximization algorithm (described by Bladt and Soerensen, 2005)

}
\usage{
gmEM(tmabs, te, gmguess, eps = 1e-06, niter = 10000, expmethod = "PadeRBS", 
verbose = FALSE)
}

\arguments{
  \item{tmabs}{
matrix of absolute transition frequencies
}
  \item{te}{
time elapsed in transition process
}
  \item{gmguess}{
initial guess (for generator matrix)
}
  \item{eps}{
stop criterion: stop, if relative change in log-likelihood is smaller than eps
}
  \item{niter}{
stop criterion: maximum number of iterations
}
  \item{expmethod}{
method for computation of matrix exponential, by default "PadeRBS" is chosen (see \code{?expm} from \code{expm} package for more information)
}
  \item{verbose}{
verbose mode
}
}
\details{
A maximum likelihood generator matrix estimate is derived by an instance of the expectation-maximization algorithm. 
}
\references{
M. Bladt and M. Soerensen: Statistical Inference for Discretely Observed Markov Jump Processes. Journal of the Royal Statistical Society B 67(3):395-410, 2005


}
\author{
Marius Pfeuffer
}
\examples{
data(tm_abs)

## Initial guess for generator matrix (absorbing default state)
gm0=matrix(1,8,8)
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0

## Derive expectation-maximization algorithm generator matrix estimate
gmem=gmEM(tmabs=tm_abs,1,gmguess=gm0,verbose=TRUE)
gmem
}

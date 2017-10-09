\name{expmMC}
\alias{expmMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix Exponential Function
}
\description{
Computation of the matrix exponential and interface to the markovchain package
}
\usage{
expmMC(gm, t, method = "PadeRBS", order = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gm}{
generator matrix (either gm or matrix object)
}
  \item{t}{
time horizon
}
  \item{method}{
method to compute the matrix exponential, see expm for details
}
  \item{order}{
order for Pade or Taylor method, see expm for details
}
}
\details{
An interface to the markovchain package is provided so that the resulting transition matrices are returned as markovchain objects and can be further processed in the markovchain package, e.g., visualized (see example below).
}

\references{
G. A. Spedicato: Discrete Time Markov Chains with R. The R Journal (To Appear), 2017
}
\author{
G. A. Spedicato, M. Pfeuffer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{expm}}
}
\examples{
### Exemplary Transition Matrix
tm0=matrix(1:16,4,4)
tm0=tm0/rowSums(tm0)
tm0[4,]=c(0,0,0,1)

### Generator Matrix Estimate
gm_est=gm(tm0,te=1,method="DA")
gm_est

### Matrix Exponential and Conversion to markovchain object 
### (markovchain package needs to be installed)
mc=expmMC(gm_est,.5)
if(require("markovchain")==TRUE){
  plot(mc)
}
}


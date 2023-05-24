\name{rNijTRiT_Unif}
\alias{rNijTRiT_Unif}
\title{
C++ Based Uniformization Sampling
}
\description{
Function for generating initial and endpoint-conditioned Markov process sampling paths for a given discrete-time transition matrix
}
\usage{
rNijTRiT_Unif(tmabs, te, gm, tpm) 
}

\arguments{
  \item{tmabs}{
matrix of absolute transition frequencies
}
  \item{te}{
time elapsed in transition process
}
  \item{gm}{
generator matrix
}
  \item{tpm}{
discrete-time transition probability matrix, matrix exponential of gm
} 
}
\details{
Function for the simulation of paths from an endpoint-conditioned Markov process. Returns number of transitions NijT and cumulative holding times RiT.
}

\value{
endpoint-conditioned sampling path
}

\references{
J. Fintzi: R Package ECctmc, 2016.

A. Hobolth and E. A. Stone: Simulation from Endpoint-Conditioned, Continuous-Time Markov Chains on a Finite State Space, with Applications to Molecular Evolution. Annals of Applied Statistics 3(3):1204-1231, 2009
}
\author{
Jon Fintzi, Marius Pfeuffer
}
\examples{
data(tm_abs)

## Generator Matrix
gm=matrix(1,8,8)
diag(gm)=0
diag(gm)=-rowSums(gm)
gm[8,]=0

## Transition Probability Matrix
library(expm)
te=1
tpm=expm(gm*te)

rNijTRiT_Unif(tm_abs,te,gm,tpm)
}

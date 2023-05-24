\name{ctmcdlogLik}
\alias{ctmcdlogLik}
\title{
Discrete-Time Data Log-Likelihood Function
}
\description{
Function for evaluating the likelihood function of a continuous-time Markov chain given discrete-time data.
}
\usage{
ctmcdlogLik(gm, tmabs, te)
}

\arguments{
\item{gm}{
generator matrix of continuous-time Markov chain
}
  \item{tmabs}{
matrix of absolute transition frequencies
}
  \item{te}{
time elapsed in transition process
}
 
}
\value{
log-likelihood value
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

## Log-likelihood of initial guess
ctmcdlogLik(gm0,tm_abs,1)
}

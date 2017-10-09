\name{ciGS}
\alias{ciGS}
\title{
Equal Tailed Credibility Interval
}
\description{
Function to derive an equal-tailed credibility interval from Gibbs sampling draws as described by Bladt and Soerensen, 2009. 
}
\usage{
ciGS(x, alpha)
}

\arguments{
  \item{x}{
Gibbs sampling output object
}
  \item{alpha}{
significance level
}
}
\details{
Computation of credibility interval es empirical quantiles of Gibbs sampling draws.
}
\references{
M. Bladt and M. Soerensen. Efficient Estimation of Transition Rates Between Credit Ratings from Observations at Discrete Time Points. Quantitative Finance 9(2):147-160, 2009
}
\author{
Marius Pfeuffer
}

\examples{
\dontrun{
data(tm_abs)

## Prior
pr=list()
pr[[1]]=matrix(1,8,8)
pr[[1]][8,]=0

pr[[2]]=c(rep(5,7),Inf)

## Derive Gibbs sampling generator matrix estimate
gmgs=gm(tm=tm_abs,te=1,method="GS",sampl_method="Unif",prior=pr,burnin=10,niter=100)

cigs=ciGS(gmgs,alpha=.05)
cigs
}
}
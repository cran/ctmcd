\name{ctmcd-package}
\alias{ctmcd-package}
\alias{ctmcd}
\docType{package}
\title{
\packageTitle{ctmcd}
}
\description{
Functions for estimating Markov generator matrices from discrete-time observations.
}
\author{
\packageAuthor{ctmcd}

Maintainer: \packageMaintainer{ctmcd}
}
\references{
M. Pfeuffer: ctmcd: An R Package for Estimating the Parameters of a Continuous-Time Markov Chain from Discrete-Time Data. The R Journal 9(2):127-141, 2017

M. Pfeuffer. Generator Matrix Approximation Based on Discrete-Time Rating Migration Data. Master
Thesis, Ludwig Maximilian University of Munich, 2016

R. B. Israel et al.: Finding Generators for Markov Chains via Empirical Transition Matrices, with Applications to Credit Ratings. Mathematical Finance 11(2):245-265, 2001

E. Kreinin and M. Sidelnikova: Regularization Algorithms for Transition Matrices. Algo Research Quarterly 4(1):23-40, 2001

M. Bladt and M. Soerensen: Statistical Inference for Discretely Observed Markov Jump Processes. Journal of the Royal Statistical Society B 67(3):395-410, 2005

}

\keyword{ package }

\examples{
\donttest{
data(tm_abs)

## Maximum Likelihood Generator Matrix Estimate
gm0=matrix(1,8,8)
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0

gmem=gm(tm_abs,te=1,method="EM",gmguess=gm0)
plot(gmem)

## Confidence Interval
ciem=gmci(gmem,alpha=0.05)
plot(ciem)
}
}

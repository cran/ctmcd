\name{gm}
\alias{gm}

\title{
Generator Matrix Estimation
}
\description{
Generic function to estimate the parameters of a continuous Markov chain
}
\usage{
gm(tm, te, method, ...)
}

\arguments{
  \item{tm}{
matrix of either absolute transition frequencies (if method is "EM" or "GS") or relative transition frequencies (if method is "DA", "WA" of "QO")
}
  \item{te}{
time elapsed in transition process
}
  \item{method}{
method to derive generator matrix: "DA" - Diagonal Adjustment, "WA" - Weighted Adjustment, "QO" - Quasi-Optimization, "EM" - Expectation-Maximization Algorithm, "GS" - Gibbs Sampler
}
 \item{...}{Additional Arguments:
 \itemize{
  \item{gmguess:}{
initial guess for generator matrix estimation procedure (if method is "EM")
}
  \item{prior:}{
prior parametrization (if method is "GS")
}
  \item{burnin:}{
burn-in period (if method is "GS")
}
  \item{eps:}{
convergence criterion (if method is "EM")
}
  \item{conv_pvalue,conv_freq:}{
convergence criterion (if method is "GS")
  }
  \item{niter:}{
maximum number of iterations (if method is "EM" or "GS")
}
  \item{sampl_func:}{
optional self-written path sampling function for endpoint-conditioned Markov processes (if method is "GS")
}
  \item{combmat:}{
matrix stating combined use of modified rejection sampling / uniformization sampling algorithms (if method is "GS")
}
  \item{sampl_method:}{
sampling method for deriving endpoint-conditioned Markov process path: "Unif" - Uniformization Sampling, "ModRej" - Modified Rejection Sampling (if method is "GS")
}
  \item{logmethod:}{
method to compute matrix logarithm (if method is "DA", "WA" or "QO", see \code{?logm} from \code{expm} package for more information)
}
  \item{expmethod:}{
method to compute matrix exponential (if method is "EM" or "GS", see \code{?expm} from \code{expm} package for more information)
}
  \item{verbose:}{
verbose mode (if method is "EM" or "GS")
}
}
}
}
\details{
The methods "DA", "WA" and "QO" provide adjustments of a matrix logarithm based candidate solution, "EM" gives the maximum likelihood estimate and "GS" a posterior mean estimate in a Bayesian setting with conjugate Gamma priors.
}


\references{
M. Pfeuffer: Generator Matrix Approximation Based on Discrete Time Rating Migration Data. Master Thesis, University of Munich, 2016

Y. Inamura: Estimating Continuous Time Transition Matrices from Discretely Observed Data. Bank of Japan Working Paper Series, 2006

R. B. Israel et al.: Finding Generators for Markov Chains via Empirical Transition Matrices, with Applications to Credit Ratings. Mathematical Finance 11(2):245-265, 2001

E. Kreinin and M. Sidelnikova: Regularization Algorithms for Transition Matrices. Algo Research Quarterly 4(1):23-40, 2001

M. Bladt and M. Soerensen: Statistical Inference for Discretely Observed Markov Jump Processes. Journal of the Royal Statistical Society B 67(3):395-410, 2005

}
\author{
Marius Pfeuffer
}


\seealso{
\code{\link{gmDA}}, \code{\link{gmWA}}, \code{\link{gmQO}}, \code{\link{gmEM}}, \code{\link{gmGS}}
}

\examples{
data(tm_abs)

## Maximum Likelihood Generator Matrix Estimate
gm0=matrix(1,8,8)
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0

gmem=gm(tm_abs,te=1,method="EM",gmguess=gm0)
gmem

## Quasi Optimization Estimate
tm_rel=rbind((tm_abs/rowSums(tm_abs))[1:7,],c(rep(0,7),1))

gmqo=gm(tm_rel,te=1,method="QO")
gmqo
}

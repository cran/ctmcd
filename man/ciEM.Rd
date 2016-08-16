\name{ciEM}
\alias{ciEM}
\title{
Wald Confidence Interval
}
\description{
Function to derive a Wald confidence interval for expectation-maximization based maximum likelihood generator matrix estimate. The Fisher information estimation method of Oakes, 1999 is employed.
}
\usage{
ciEM(x, alpha, eps = 1e-04, expmethod = "PadeRBS")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
expectation-maximization output object
}
  \item{alpha}{
significance level
}
  \item{eps}{
threshold for which generator matrix parameters are assumed to be fixed at zero
}
  \item{expmethod}{
method for computation of matrix exponential, by default "PadeRBS" is chosen (see \code{?expm} from \code{expm} package for more information)
}
}
\details{
Computation of confidence interval on the basis of the last iteration of the expectation-maximization algorithm.
}
\references{
M. Bladt and M. Soerensen. Efficient Estimation of Transition Rates Between Credit Ratings from Observations at Discrete Time Points. Quantitative Finance 9(2):147-160, 2009

D. Oakes. Direct calculation of the information matrix via the EM algorithm. Journal of the Royal Statistical Society B 61(2):479-482, 1999
}
\author{
Marius Pfeuffer
}
\examples{
\dontrun{
data(tm_abs)
## Initial guess for generator matrix (absorbing default state)
gm0=matrix(1,8,8)
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0

## Derive expectation-maximization algorithm generator matrix estimate
gmem=gm(tm=tm_abs,te=1,method="EM",gmguess=gm0,verbose=TRUE)

## derive confidence interval
ciem=ciEM(gmem, alpha=0.05)
ciem
}
}
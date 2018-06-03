\name{tmci}
\alias{tmci}
\title{
Delta Method Confidence Intervals for Matrix Exponential Transformations of Generator Matrix Objects
}
\description{
Generic function to derive delta method based confidence intervals for matrix exponential transformations of "EM" based generator matrix objects
}
\usage{
tmci(gmem, alpha, te, eps = 1e-04, expmethod = "PadeRBS")
}
\arguments{
  \item{gmem}{
an "EM" generator matrix object
}
  \item{alpha}{
significance level
}
  \item{te}{
discrete time horizon for which the interval is supposed to be computed
}
 \item{eps}{
threshold for which generator matrix parameters are assumed to be fixed at zero
}
  \item{expmethod}{
method to compute matrix exponentials (see \code{?expm} from \code{expm} package for more information)
}
}
\details{
Confidence intervals for discrete-time transition matrix predictions given generator matrix estimates are computed by using the delta method for matrix exponential transformations. 
}

\references{
G. dos Reis, M. Pfeuffer, G. Smith: Capturing Rating Momentum in the Estimation of Probabilities of Default, With Application to Credit Rating Migrations (In Preparation), 2018
}

\examples{
\dontrun{
data(tm_abs)

## Maximum Likelihood Generator Matrix Estimate
gm0=matrix(1,8,8)
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0

gmem=gm(tm_abs,te=1,method="EM",gmguess=gm0)

## 2.5 Year Transition Matrix Confidence Interval
citm=tmci(gmem,alpha=0.05,te=2.5)
citm
}
}
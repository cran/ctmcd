\name{gmci.default}
\alias{gmci.default}
\alias{ciEMBS}
\alias{ciEMSdR}
\alias{ciEMoFI}
\alias{ciEMDM}
\alias{ciGS}
\title{
Confidence / Credibility Intervals for Generator Matrix Objects
}
\description{
Default function to derive confidence / credibility intervals for "EM" or "GS" based generator matrix objects
}
\usage{
\method{gmci}{default}(gm, alpha, eps = 1e-04, cimethod="Direct", expmethod = "PadeRBS", ...)
}
\arguments{
  \item{gm}{
a "EM" or "GS" generator matrix object
}
  \item{alpha}{
significance level
}
  \item{eps}{
threshold for which generator matrix parameters are assumed to be fixed at zero (if "EM" object)
}
 \item{cimethod}{
"Direct" or "SdR" use analytical expressions of the Fisher information matrix, "BS" emloy the numerical expressions of Bladt and Soerensen, 2009 (if "EM" object)
}
  \item{expmethod}{
method to compute matrix exponentials (see \code{?expm} from \code{expm} package for more information)
}
  \item{\dots}{
additional arguments
}
}
\details{
If gm is based on the "EM" method (expectation-maximization algorithm), the function computes a Wald confidence interval based on the method of Oakes, 1999. IF gm is based on the "GS" method (Gibbs sampler), the function computes an equal-tailed credibility interval.
}

\references{

G. dos Reis, M. Pfeuffer, G. Smith: Capturing Rating Momentum in the Estimation of Probabilities of Default, With Application to Credit Rating Migrations (In Preparation), 2018
}

\author{
Marius Pfeuffer
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

## Oakes Confidence Interval
ciem=gmci(gmem,alpha=0.05)
ciem
}
}
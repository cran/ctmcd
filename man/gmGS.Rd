\name{gmGS}
\alias{gmGS}

\title{
Gibbs Sampler
}
\description{
Function for deriving a Markov generator matrix estimate by Gibbs sampling (described by Bladt and Soerensen, 2005)
}
\usage{
gmGS(tmabs, te, prior, burnin, conv_pvalue = 0, conv_freq = 10,
niter = 10000, sampl_method = "Unif", expmethod = "PadeRBS", verbose = FALSE, 
combmat=NULL, sampl_func = NULL)
}
\arguments{
  \item{tmabs}{
matrix of absolute transition frequencies
}
  \item{te}{
time elapsed in transition process
}
  \item{prior}{
list of prior parameters (Gamma prior)
}
  \item{burnin}{
number of burn-in iterations
}
  \item{conv_pvalue}{
convergence criterion: stop, if Heidelberger and Welch's diagnostic assumes convergence (see coda package), convergence check is only employed if conv_pvalue>0
}
 \item{conv_freq}{
convergence criterion: absolute frequency of convergence evaluations
}
  \item{niter}{
stop criterion: stop, if maximum number of iterations is exceeded
}
  \item{sampl_method}{
method for sampling paths from endpoint-conditioned Markov processes. options: "Unif" - Uniformization sampling, "ModRej" - Modified Rejection Sampling
}
  \item{expmethod}{
method for computation of matrix exponential, by default "PadeRBS" is chosen (see \code{?expm} from \code{expm} package for more information)
}
  \item{verbose}{
verbose mode
}
  \item{combmat}{
matrix specifying the combined use of sampling methods: "U" - uniformization sampling, "M" - modified rejection sampling
}
  \item{sampl_func}{
interface for own endpoint-conditioned Markov process sampling function
}
}
\details{
A posterior mean generator matrix estimate is derived by Gibbs Sampling. The gamma distribution is used as prior.
}
 
\references{
M. Bladt and M. Soerensen: Statistical Inference for Discretely Observed Markov Jump Processes. Journal of the Royal Statistical Society B 67(3):395-410, 2005


}
\author{
Marius Pfeuffer
}

\seealso{
\code{\link{rNijTRiT_ModRej}}, \code{\link{rNijTRiT_Unif}}
}

\examples{
data(tm_abs)

## Example prior parametrization (absorbing default state)
pr=list()
pr[[1]]=matrix(1,8,8)
pr[[1]][8,]=0

pr[[2]]=c(rep(5,7),Inf)

## Derive Gibbs sampling generator matrix estimate
\dontrun{
gmgs=gmGS(tmabs=tm_abs,te=1,sampl_method="Unif",prior=pr,burnin=10,niter=100,verbose=TRUE)
gmgs
}

}

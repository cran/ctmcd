\name{gmQO}
\alias{gmQO}
\title{
Quasi-Optimization
}
\description{
Function for deriving a Markov generator matrix estimate based on the quasi-optimization procedure of Kreinin and Sidelnikova, 2001
}
\usage{
gmQO(tmrel, te, logmethod = "Eigen")
}
\arguments{
  \item{tmrel}{
matrix of relative transition frequencies
}
  \item{te}{
time elapsed in transition process
}
  \item{logmethod}{
method for computation of matrix logarithm, by default eigendecomposition is chosen (see \code{?logm} from \code{expm} package for more information)
}
}
\details{
From the set of possible Markov generator matrices, the one is chosen which is closest to a matrix logarithm based candidate solution in terms of sum of squared deviations.
}

\value{
generator matrix estimate
}
\references{
E. Kreinin and M. Sidelnikova: Regularization Algorithms for Transition Matrices. Algo Research Quarterly 4(1):23-40, 2001
}
\author{
Marius Pfeuffer
}
\examples{
data(tm_abs)
## Derive matrix of relative transition frequencies
data(tm_abs)
tm_rel=rbind((tm_abs/rowSums(tm_abs))[1:7,],c(rep(0,7),1))

## Derive quasi optimization generator matrix estimate
gmqo=gmQO(tm_rel,1)
gmqo
}
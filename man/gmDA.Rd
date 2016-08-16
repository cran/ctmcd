\name{gmDA}
\alias{gmDA}
\title{
Diagonal Adjustment
}
\description{
Function for deriving a Markov generator matrix estimate based on the diagonal adjustment method of Israel et al., 2001
}
\usage{
gmDA(tmrel, te, logmethod = "Eigen")
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
A candidate solution is derived by the matrix logarithm and then adjusted in order to fulfil the properties of a Markov generator matrix.
}
\references{
R. B. Israel et al.: Finding Generators for Markov Chains via Empirical Transition Matrices, with Applications to Credit Ratings. Mathematical Finance 11(2):245-265, 2001
}

\author{
Marius Pfeuffer
}

\examples{
## Derive matrix of relative transition frequencies
data(tm_abs)
tm_rel=rbind((tm_abs/rowSums(tm_abs))[1:7,],c(rep(0,7),1))

## Derive diagonal adjustment generator matrix estimate
gmda=gmDA(tm_rel,1)
gmda
}

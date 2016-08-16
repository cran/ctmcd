\name{tm_abs}
\alias{tm_abs}
\docType{data}
\title{
Single Year Corporate Credit Rating Transititions
}
\description{
Matrix of Standard and Poor's Global Corporate Rating Transition Frequencies 2000 (NR Removed)
}
\usage{data("tm_abs")}
\format{
  The format is:
 num [1:8, 1:8] 17 2 0 0 0 0 0 0 1 455 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:8] "AAA" "AA" "A" "BBB" ...
  ..$ : chr [1:8] "AAA" "AA" "A" "BBB" ...
}

\references{
European Securities and Markets Authority, 2016 

https://cerep.esma.europa.eu/cerep-web/statistics/transitionMatrice.xhtml
}
\examples{
data(tm_abs)

## Matrix of relative transition frequencies
tm_rel=rbind((tm_abs/rowSums(tm_abs))[1:7,],c(rep(0,7),1))
tm_rel

}
\keyword{datasets}

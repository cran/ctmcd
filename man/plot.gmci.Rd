\name{plot.gmci}
\alias{plot.gmci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Function for Generator Matrix Confidence / Credibility Interval Objects
}
\description{
Function for visualizing the boundaries of generator matrix confidence / credibility intervals
}
\usage{
\method{plot}{gmci}(x, mattext, col = c("grey", "red"), main, las = 1, xlab = "To",
ylab = "From", xnames, ynames, cex = 1, fig = 2, opacity_factor, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a generator matrix confidence / credibility interval object
}
  \item{mattext}{
optional: matrix of strings replacing the parameter estimates
}
  \item{col}{
two element vector of basis colors for positive and negative parameter estimate entries
}
  \item{main}{
optional: plot title
}
  \item{las}{
orientation of x and y axis elements
}
  \item{xlab}{
x axis name
}
  \item{ylab}{
y axis name
}
  \item{xnames}{
description of x axis elements
}
  \item{ynames}{
description of y axis elements
}
  \item{cex}{
font size
}
  \item{fig}{
number of significant figures to be plotted
}
  \item{opacity_factor}{
two element vector for specification of opacity for positive and negative parameter entry highlighting (must be greater than zero)

}
  \item{\dots}{
additional arguments
}
}
\author{
Marius Pfeuffer
}
\seealso{
\code{\link{print.gmci}}, \code{\link{plotM}}
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
plot(gmem)

## Confidence Interval
ciem=gmci(gmem,alpha=0.05)
plot(ciem)
}
}

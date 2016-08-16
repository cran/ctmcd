\name{plot.gm}
\alias{plot.gm}

\title{
Plot Function for Generator Matrix Estimation Objects
}
\description{
Function for visualizing the output of a generator matrix estimation procedure.
}
\usage{
\method{plot}{gm}(x, mattext, col = c("grey", "red"), main = x$method, las = 1, 
xlab = "To", ylab = "From", xnames, ynames, cex = 1, fig = 3, opacity_factor, ...)
}
\arguments{
  \item{x}{
a generator matrix estimation object
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
number of significant figure to be plotted
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
\code{\link{print.gm}}, \code{\link{summary.gm}}, \code{\link{plotM}}
}
\examples{
data(tm_abs)

## Maximum Likelihood Generator Matrix Estimate
gm0=matrix(1,8,8)
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0

gmem=gm(tm_abs,te=1,method="EM",gmguess=gm0)
plot(gmem)
}

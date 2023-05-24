\name{plotM}
\alias{plotM}
\title{
Matrix Plot Function
}
\description{
Function to visualize matrices
}
\usage{
plotM(mat, mattext, col = c("grey", "red"), main, las = 1, xlab = "To", 
ylab = "From", xnames, ynames, cex = min(1, nrow(mat)/8), fig = 3, opacity_factor)
}
\arguments{
  \item{mat}{
a matrix
}
  \item{mattext}{
optional: matrix of strings replacing the original matrix entries
}
  \item{col}{
two element vector of basis colors for positive and negative matrix entries
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
}

\value{
no value, plot function
}

\author{
Marius Pfeuffer
}\seealso{
\code{\link{plot.gm}}, \code{\link{plot.gmci}}
}
\examples{gm0=matrix(1,8,8)
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0

plotM(gm0)
}
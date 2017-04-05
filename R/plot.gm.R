plot.gm <-
function (x, mattext, col = c("grey", "red"), main = x$method, 
    las = 1, xlab = "To", ylab = "From", xnames, ynames, cex = 1, 
    fig = 3, opacity_factor, ...) 
{
    plotM(x$par, mattext, col, main, las, xlab, ylab, xnames, 
        ynames, cex, fig, opacity_factor)
}

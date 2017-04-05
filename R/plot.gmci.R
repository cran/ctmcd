plot.gmci <-
function (x, mattext, col = c("grey", "red"), main, las = 1, 
    xlab = "To", ylab = "From", xnames, ynames, cex = 1, fig = 2, 
    opacity_factor, ...) 
{
    if (missing(main)) {
        main = paste0(round(100 * (1 - x$alpha), 3), "% ", x$method)
    }
    n = nrow(x$par)
    if (missing(mattext)) {
        mattext = matrix(paste0("[", signif(x$lower, fig), "; ", 
            signif(x$upper, fig), "]"), n, n)
        mattext[which(is.na(x$lower) & is.na(x$upper))] = ""
    }
    plotM(x$par, mattext, col, main, las, xlab, ylab, xnames, 
        ynames, cex, fig, opacity_factor)
}

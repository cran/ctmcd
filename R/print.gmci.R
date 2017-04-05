print.gmci <-
function (x, ...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\nMethod:\n")
    print(x$method)
    cat("\nLower Limit:\n")
    print(signif(x$lower, 4))
    cat("\nUpper Limit:\n")
    print(signif(x$upper, 4))
}

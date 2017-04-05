ciGS <-
function (x, alpha) 
{
    draws = length(x$draws)
    n = nrow(x$draws[[1]])
    draws = length(x$draws)
    data = matrix(unlist(x$draws), n^2, draws)
    lowermat = matrix(0, n, n)
    uppermat = matrix(0, n, n)
    for (row in 1:n) {
        for (col in 1:n) {
            lowermat[row, col] = quantile(data[(col - 1) * n + 
                row, ], alpha/2)
            uppermat[row, col] = quantile(data[(col - 1) * n + 
                row, ], 1 - alpha/2)
        }
    }
    limits = list(lower = lowermat, upper = uppermat)
    limits
}

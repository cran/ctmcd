gmWA <-
function (tmrel, te, logmethod = "Eigen") 
{
    n = nrow(tmrel)
    gmest = logm(tmrel, method = logmethod)/te
    for (i in 1:(n - 1)) {
        gmiNeg = -sum((gmest * (gmest < 0))[i, setdiff(1:n, i)])
        gmiPos = sum((gmest * (gmest > 0))[i, ])
        gmest[i, setdiff(1:n, i)] = gmest[i, setdiff(1:n, i)] - 
            gmiNeg/gmiPos * abs(gmest[i, setdiff(1:n, i)])
    }
    gmest[setdiff(which(gmest < 0), seq(1, n^2, n + 1))] = 0
    return(gmest)
}

gmDA <-
function (tmrel, te, logmethod = "Eigen") 
{
    gmest = logm(tmrel, method = logmethod)/te
    gmest[which(gmest < 0)] = 0
    diag(gmest) = -rowSums(gmest)
    return(gmest)
}

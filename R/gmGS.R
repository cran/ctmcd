gmGS <-
function (tmabs, te, prior, burnin, conv_pvalue = 0, conv_freq = 10, 
    niter = 10000, sampl_method = "Unif", expmethod = "PadeRBS", 
    verbose = FALSE, combmat = NULL, sampl_func = NULL) 
{
    n = ncol(tmabs)
    gmest = matrix(0, n, n)
    for (i in 1:n) {
        for (j in setdiff(1:n, i)) {
            gmest[i, j] = rgamma(1, prior[[1]][i, j], rate = prior[[2]][i])
        }
    }
    diag(gmest) = -rowSums(gmest)
    gmestList = list()
    count = 0
    if(verbose==TRUE){
      pb=txtProgressBar(min=0,max=burnin+niter,title="Burn-in",style=3)
    }
    repeat {
        if (sampl_method == "Unif") {
            draw = rNijTRiT_Unif(tmabs, te, gmest, expm(gmest, 
                method = expmethod))
        }
        else if (sampl_method == "ModRej") {
            draw = rNijTRiT_ModRej(tmabs, te, gmest)
        }
        else if (sampl_method == "Comb") {
            Mmat = (combmat == "U") * tmabs
            Umat = (combmat == "M") * tmabs
            Mdraw = rNijTRiT_ModRej(Mmat, te, gmest)
            Udraw = rNijTRiT_Unif(Umat, te, gmest, expm(gmest, 
                method = expmethod))
            draw = list()
            draw$NijT = Mdraw$NijT + Udraw$NijT
            draw$RiT = Mdraw$RiT + Udraw$RiT
        }
        else if (sampl_method == "Ext") {
            draw = sampl_func(tmabs, te, gmest)
        }
        for (i in 1:n) {
            for (j in setdiff(1:n, i)) {
                gmest[i, j] = rgamma(1, draw$NijT[i, j] + prior[[1]][i, 
                  j], rate = (draw$RiT[i] + prior[[2]][i]))
            }
        }
        diag(gmest) = 0
        diag(gmest) = -rowSums(gmest)
        count = count + 1
        if (count > burnin) {
            gmestList[[length(gmestList) + 1]] = gmest
        }
        if (verbose == TRUE) setTxtProgressBar(pb,count)
        if (conv_pvalue>0 && count > burnin + 10 && (count - burnin)%%round(niter/conv_freq) == 
            0) {
            parmat = t(matrix(unlist(gmestList), n^2))
            parmat2 = parmat[, which(colSums(parmat) > 0)]
            mcmcobj = as.mcmc(parmat2)
            heidel = heidel.diag(mcmcobj, pvalue = conv_pvalue)
            if (verbose == TRUE) {
                print("Check Convergence")
            }
            if (sum(heidel[, 1]) == length(heidel[, 1]) & sum(heidel[, 
                4]) == length(heidel[, 4])) {
                if (verbose == TRUE) {
                  print("Converged according to Heidelberger and Welch Criterion")
                }
                break
            }
            if (verbose == TRUE) {
                print("Not Converged according to Heidelberger and Welch Criterion")
            }
        }
        if (count == (burnin + niter)) break
    }
    gmall = Reduce("+", gmestList)/(count - burnin)
    colnames(gmall) = colnames(tmabs)
    rownames(gmall) = rownames(tmabs)
    est = list(par = gmall, burnin = burnin, niter = count, draws = gmestList)
    return(est)
}

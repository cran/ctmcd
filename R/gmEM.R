gmEM <-
function (tmabs, te, gmguess, eps = 1e-06, niter = 10000, expmethod = "PadeRBS", 
                verbose = FALSE)
  
  ### Implementation according to Bladt and Sorensen
  ### Statistical Inference for Discretely Observed Markov Processes
  ### Journal of the Royal Statistical Society B, 2005
  ### and Inamura
  ### Estimating Continuous Time Transition Matrices From Discretely Observed Data
  ### Bank of Japan Working Paper Series 2006
  
{
  n = ncol(tmabs)
  if (nrow(tmabs) == n - 1) {
    tmabs = rbind(tmabs, rep(0, n))
  }
  llvec = NULL
  gmprev = gmguess
  for (k in 1:niter) {
    dtt = expm(gmprev * te, method = expmethod)
    ERiT = rep(0, n)
    for (i in 1:n) {
      uitui = matrix(0, n, n)
      uitui[i, i] = 1
      augmat = rbind(cbind(gmprev, uitui), cbind(matrix(0, 
                                                        n, n), gmprev))
      caseRiT = expm(te * augmat, method = expmethod)[1:n, 
                                                      (n + 1):(2 * n)]
      ERiT[i] = ERiT[i] + sum(tmabs * caseRiT/dtt, na.rm = TRUE)
    }
    ENijT = matrix(0, n, n)
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {
        uituj = matrix(0, n, n)
        uituj[i, j] = 1
        augmat = rbind(cbind(gmprev, uituj), cbind(matrix(0, 
                                                          n, n), gmprev))
        caseNijT = gmprev[i, j] * expm(te * augmat, method = expmethod)[1:n, 
                                                                        (n + 1):(2 * n)]
        ENijT[i, j] = ENijT[i, j] + sum(tmabs * caseNijT/dtt, 
                                        na.rm = TRUE)
      }
    }
    gmest = rbind(ENijT/ERiT)
    diag(gmest) = -rowSums(gmest)
    gmprev = gmest
    ll=ctmcdlogLik(gmest,tmabs,te)
    llvec = c(llvec, ll)
    rel_err = abs(llvec[length(llvec)] - llvec[length(llvec) - 
                                                 1])/abs(llvec[length(llvec) - 1])
    if (verbose == TRUE) {
      print(paste0("Iteration: ", k, "  Log-Likelihood: ", 
                   signif(ll, 6)))
    }
    if (k > 2 && rel_err <= eps) {
      break
    }
  }
  colnames(gmest) = colnames(tmabs)
  rownames(gmest) = rownames(tmabs)
  est = list(par = gmest, niter = k, eps = rel_err, ll = llvec, 
             ENijT = ENijT, ERiT = ERiT)
  return(est)
}

gmQO <-
function(tmrel, te, logmethod = "Eigen") {
  ### Implementation according to Kreinin and Sidelnikova
  ### Regularization Algorithms for Transition Matrices
  ### Algo Research Quarterly 4(1):23-40, 2001
  n = nrow(tmrel)
  gmest = logm(tmrel, method = logmethod) / te
  for (i in 1:n) {
    a = gmest[i,]
    lambda = mean(a)
    aorder = order(a + lambda)
    aest = a[aorder]
    for (m in 2:n) {
      if ((n - m + 1) * aest[m + 1] - (aest[1] + sum(aest[(m + 1):n])) >= 0) {
        mstar = m
        break
      }
    }
    zstar = NULL
    for (j in 1:n) {
      if (j %in% 2:mstar) {
        zstar = c(zstar, 0)
      } else{
        zstar = c(zstar, aest[j] - 1 / (n - mstar + 1) * (aest[1] + sum(aest[(mstar +
                                                                                1):n])))
      }
    }
    gmest[i,] = zstar[order(aorder)]
  }
  return(gmest)
}

gmDA <-
function(tmrel, te, logmethod = "Eigen") {
  ### Implementation accoring to Israel et al.
  ### Finding Generators for Markov Chains via Empirical Transition Matrices,
  ### With Applications to Credit Ratings
  ### Mathematical Finance 11(2):245-265, 2001
  gmest = logm(tmrel, method = logmethod) / te
  gmest[which(gmest < 0)] = 0
  diag(gmest) = -rowSums(gmest)
  return(gmest)
}

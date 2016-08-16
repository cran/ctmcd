gmGS <-
function(tmabs,
                te,
                prior,
                burnin,
                eps = 1e-6,
                niter = 1e4,
                sampl_method = "Unif",
                expmethod = "PadeRBS",
                verbose = FALSE,
                combmat = NULL,
                sampl_func = NULL) {
  ### Implementation according to Bladt and Sorensen
  ### Statistical Inference for Discretely Observed Markov Processes
  ### Journal of the Royal Statistical Society B, 2005
  ### and Inamura
  ### Estimating Continuous Time Transition Matrices From Discretely Observed Data
  ### Bank of Japan Working Paper Series 2006
  
  n = ncol(tmabs)
  
  ### Initialization: Draw from Prior
  gmest = matrix(0, n, n)
  for (i in 1:n) {
    for (j in setdiff(1:n, i)) {
      gmest[i, j] = rgamma(1, prior[[1]][i, j], rate = prior[[2]][i])
    }
  }
  diag(gmest) = -rowSums(gmest)
  
  l1metric = NULL
  gmall = matrix(0, n, n)
  gmestList = list()
  count = 0
  
  repeat {
    ### Draw NijT and RiT From Missing Data Full Conditional: Simulate Paths
    if (sampl_method == "Unif") {
      draw = rNijTRiT_Unif(tmabs, te, gmest, expm(gmest, method = expmethod))
    } else if (sampl_method == "ModRej") {
      draw = rNijTRiT_ModRej(tmabs, te, gmest)
    } else if (sampl_method == "Comb") {
      Mmat = (combmat == "U") * tmabs
      Umat = (combmat == "M") * tmabs
      Mdraw = rNijTRiT_ModRej(Mmat, te, gmest)
      Udraw = rNijTRiT_Unif(Umat, te, gmest, expm(gmest, method = expmethod))
      draw = list()
      draw$NijT = Mdraw$NijT + Udraw$NijT
      draw$RiT = Mdraw$RiT + Udraw$RiT
    } else if (sampl_method == "Ext") {
      draw = sampl_func(tmabs, te, gmest)
    }
    
    ### Draw Generator Matrix from Parameter Full Conditional
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {
        gmest[i, j] = rgamma(1, draw$NijT[i, j] + prior[[1]][i, j], rate = (draw$RiT[i] +
                                                                              prior[[2]][i]))
      }
    }
    diag(gmest) = 0
    diag(gmest) = -rowSums(gmest)
    
    count = count + 1
    
    if (count > burnin) {
      gmestList[[length(gmestList) + 1]] = gmest
      if (count > burnin + 1) {
        gmall_prev = gmall
      }
      gmall = gmall + gmest
      if (count > burnin + 1) {
        l1metric = c(l1metric, sum(abs(
          gmall / (count - burnin) - gmall_prev / (count - burnin - 1)
        )) / n ^ 2)
      }
    }
    
    if (verbose == TRUE) {
      if (count <= burnin) {
        print(paste0("Burn-in Iteration: ", count))
      } else if (count == burnin + 1) {
        print(paste0("Iteration: ", count - burnin))
      } else if (count > burnin + 1) {
        print(paste0(
          "Iteration: ",
          count - burnin,
          "  L1 Metric: ",
          signif(l1metric[length(l1metric)], 6)
        ))
      }
    }
    
    if (count > burnin + 2) {
      if (sum(l1metric[c(count - burnin - 2, count - burnin - 1)] < eps) == 2) {
        break
      }
    }
    
    if (count == (burnin + niter)) {
      break
    }
    
  }
  colnames(gmall) = colnames(tmabs)
  rownames(gmall) = rownames(tmabs)
  est = list(
    par = gmall / (count - burnin),
    burnin = burnin,
    niter = count,
    l1metric = l1metric,
    draws = gmestList
  )
  return(est)
}

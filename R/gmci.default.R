gmci.default <-
  function (gm, alpha, eps = 1e-04, cimethod= "Direct", expmethod = "PadeRBS", ...) 
  {
    if (gm$method == "Expectation-Maximization Algorithm") {
      if(cimethod=="SdR"){
        limits = ciEMSdR(gm, alpha, eps, expmethod)
      }
      if(cimethod=="BS"){
        limits = ciEMBS(gm, alpha, eps, expmethod)
      }
      if(cimethod=="Direct"){
        limits = ciEMoFI(gm, alpha, eps, expmethod)
      }
      if(cimethod=="DM"){
        limits = ciEMDM(gm, alpha, eps, expmethod)
      }
      limits$method = "Wald Confidence Interval (Oakes Standard Error)"
    }
    if (gm$method == "Gibbs Sampler") {
      limits = ciGS(gm, alpha)
      limits$method = "Equal Tailed Credibility Interval"
    }
    limits$par = gm$par
    limits$alpha = alpha
    limits$call = match.call()
    class(limits) = "gmci"
    limits
  }

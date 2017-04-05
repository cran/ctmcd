gm.default <-
function (tm, te, method, gmguess = NULL, prior = NULL, burnin = NULL, 
          eps = 1e-06, conv_pvalue=0.05, conv_freq = 10, niter = 10000, sampl_func = NULL, combmat = NULL, 
          sampl_method = "Unif", logmethod = "Eigen", expmethod = "PadeRBS", 
          verbose = FALSE, ...) 
{
  tm = as.matrix(tm)
  te = as.numeric(te)
  est = list()
  if (method == "DA") {
    est$par = gmDA(tmrel = tm, te, logmethod = logmethod)
    est$method = "Diagonal Adjustment"
  }
  else if (method == "WA") {
    est$par = gmWA(tmrel = tm, te, logmethod = logmethod)
    est$method = "Weighted Adjustment"
  }
  else if (method == "QO") {
    est$par = gmQO(tmrel = tm, te, logmethod = logmethod)
    est$method = "Quasi Optimization"
  }
  else if (method == "EM") {
    est = gmEM(tmabs = tm, te, gmguess, eps, niter, expmethod, 
               verbose)
    est$method = "Expectation-Maximization Algorithm"
  }
  else if (method == "GS") {
    est = gmGS(tmabs = tm, te, prior, burnin, conv_pvalue, conv_freq, niter, 
               sampl_method, expmethod, verbose, combmat, sampl_func)
    est$method = "Gibbs Sampler"
  }
  est$call = match.call()
  est$tm = tm
  est$te = te
  class(est) = "gm"
  est
}

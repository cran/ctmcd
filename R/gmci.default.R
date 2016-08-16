gmci.default <-
function(gm,
                        alpha,
                        eps = 1e-4,
                        expmethod = "PadeRBS",
                        ...) {
  if (gm$method == "Expectation-Maximization Algorithm") {
    limits = ciEM(gm, alpha, eps, expmethod)
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

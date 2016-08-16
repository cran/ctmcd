summary.gm <-
function(object, ...) {
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:\n")
  print(object$method)
  cat("\nParameters:\n")
  print(signif(object$par, 4))
  if (!is.null(object$burnin)) {
    cat("\nBurn-in:\n")
    print(object$burnin)
  }
  if (!is.null(object$niter)) {
    cat("\nNumber of Iterations:\n")
    print(object$niter)
  }
}

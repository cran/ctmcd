print.gm <-
function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nParameters:\n")
  print(signif(x$par, 4))
}

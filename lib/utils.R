# utils.R
# Various utility functions

#' Numerical gradient function for testing analytical
numgrad <- function(fcn, p, delta = 1e-10) {
  val1 <- fcn(p)
  
  toadd <- diag(rep(delta, length(p)))
  newp <- lapply(1:length(p), function(x) p + toadd[, x])
  newvals <- vapply(newp, fcn, numeric(1))
  grad <- (newvals - val1) / delta
  grad
}

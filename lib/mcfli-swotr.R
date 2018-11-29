# functions ported from mcfli-swotr


# mcfli-swotr/lib/beta.R --------------------------------------------------

# Symmetric finite difference functions that preserve dimension
# from 0417
findif_x <- function(dawgmat) {
  difmat <- apply(dawgmat, 2, diff)
  out1 <- rbind(difmat[1, ], difmat)
  out2 <- rbind(difmat, difmat[nrow(difmat), ])
  out <- (out1 + out2) / 2
  out
}

findif_t <- function(dawgmat) {
  difmat <- t(apply(dawgmat, 1, diff))
  if (ncol(dawgmat) == 2) {
    difmat <- t(difmat)
  }
  out1 <- cbind(difmat[, 1], difmat)
  out2 <- cbind(difmat, difmat[, ncol(difmat)])
  out <- (out1 + out2) / 2
  out
}

boot.sxenv <- function (X, Y, u, R, B) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  fit <- sxenv(X, Y, u, R, asy = F)
  muhat <- fit$muY - crossprod(fit$beta, fit$muX)
  Yfit <- matrix(1, n, 1) %*% t(muhat) + X %*% fit$beta
  res <- Y - Yfit
  bootenv <- function(i) {
    res.boot <- res[sample(1:n, n, replace = T), ]
    Y.boot <- Yfit + res.boot
    return(c(sxenv(X, Y.boot, u, R, asy = F)$beta))
  }
  bootbeta <- lapply(1:B, function(i) bootenv(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)
  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = p)
  return(bootse)
}
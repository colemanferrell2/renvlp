boot.senv <- function (X, Y, u, B) 
{
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  fit <- senv(X, Y, u, asy = F)
  Yfit <- matrix(1, n, 1) %*% t(fit$mu) + X %*% t(fit$beta)
  res <- Y - Yfit
  bootenv <- function(i) {
    res.boot <- res[sample(1:n, n, replace = T), ]
    Y.boot <- Yfit + res.boot
    return(c(senv(X, Y.boot, u, asy = F)$beta))
  }
  bootbeta <- lapply(1:B, function(i) bootenv(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)
  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = r)
  return(bootse)
}
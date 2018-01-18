boot.stenv <- function(X, Y, q, u, B){
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  fit <- stenv(X, Y, q, u, asy = F)
  Yfit <- matrix(1, n, 1) %*% t(fit$mu) + X %*% fit$beta
  res <- Y - Yfit
  bootstenv <- function(i) {
    res.boot <- res[sample(1:n, n, replace = T), ]
    Y.boot <- Yfit + res.boot
    return(c(stenv(X, Y.boot, q, u, asy = F)$beta))
  }
  bootbeta <- lapply(1:B, function(i) bootstenv(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)
  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = p)
  return(bootse)
}
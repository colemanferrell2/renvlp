boot.penv <- function(X1, X2, Y, u, B) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  
  fit <- penv(X1, X2, Y, u, asy = F)
  Yfit <- matrix(1, n, 1) %*% t(fit$mu) + X1 %*% t(fit$beta1) + X2 %*% t(fit$beta2)
  res <- Y - Yfit
  
  bootenv <- function(i) {
    res.boot <- res[sample(1:n, n, replace = T), ]
    Y.boot <- Yfit + res.boot
    return(c(penv(X1, X2, Y.boot, u, asy = F)$beta1))
  }
  
  bootbeta <- lapply(1:B, function(i) bootenv(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)

  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = r)
  return(bootse)
  
}

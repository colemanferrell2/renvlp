boot.rrenv <- function(X, Y, u, d, B) {
  
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  fit <- rrenv(X, Y, u, d, asy = F)
  Yfit <- matrix(1, n, 1) %*% t(fit$mu) + X %*% t(fit$beta)
  res <- Y - Yfit
  
  bootenvapw <- function(i) {
    res.boot <- res[sample(1:n, n, replace = T), ]
    Y.boot <- Yfit + res.boot
    return(c(rrenv(X, Y.boot, u, d, asy = F)$beta))
  }
  
  bootbeta <- lapply(1:B, function(i) bootenvapw(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)

  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = r)
  return(bootse)
  
}



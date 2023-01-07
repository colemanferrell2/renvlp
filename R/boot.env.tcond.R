boot.env.tcond <- function(X, Y, u, df, B) {
  
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  fit <- env.tcond(X, Y, u, df, asy = F)
  Yfit <- matrix(1, n, 1) %*% t(fit$mu) + X %*% t(fit$beta)
  res <- Y - Yfit
  
  bootenvtcond <- function(i) {
    res.boot <- res[sample(1:n, n, replace = T), ]
    Y.boot <- Yfit + res.boot
    return(c(env.tcond(X, Y.boot, u, df, asy = F)$beta))
  }
  
  bootbeta <- lapply(1:B, function(i) bootenvtcond(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)

  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = r)
  return(bootse)
  
}



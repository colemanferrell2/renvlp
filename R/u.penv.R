u.penv <- function(X1, X2, Y, alpha = 0.01) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  
  loglik.seq <- unlist(lapply(0:r, function(x) penv(X1, X2, Y, x, asy = F)$loglik))
  npara.seq <- r + r * p2 + r * (r + 1) / 2 + p1 * (0:r) 
  
  aic.seq <- -2 * loglik.seq + 2 * npara.seq
  bic.seq <- -2 * loglik.seq + log(n) * npara.seq
  
  u.aic <- which.min(aic.seq) - 1
  u.bic <- which.min(bic.seq) - 1
  
  lrt.test <- stats::pchisq(2 * (loglik.seq[r + 1] - loglik.seq[1:r]), npara.seq[r + 1] - npara.seq[1:r], lower.tail = F)
  
  if (any(lrt.test > alpha)) {
    u.lrt <- which(lrt.test > alpha)[1] - 1
  } else {
    u.lrt <- r
  }
  
  return(list(u.aic = u.aic, u.bic = u.bic, u.lrt = u.lrt, loglik.seq = loglik.seq, aic.seq = aic.seq, bic.seq = bic.seq))
}

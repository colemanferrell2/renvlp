u.xenv <- function(X, Y, alpha = 0.01) {
  
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  loglik.seq <- unlist(lapply(0:p, function(x) xenv(X, Y, x, asy = F)$loglik))
  npara.seq <- r + r * (r + 1) / 2 + p * (p + 1) / 2 + r * (0:p)
  
  aic.seq <- -2 * loglik.seq + 2 * npara.seq
  bic.seq <- -2 * loglik.seq + log(n) * npara.seq
  
  u.aic <- which.min(aic.seq) - 1
  u.bic <- which.min(bic.seq) - 1
  
  lrt.test <- stats::pchisq(2 * (loglik.seq[p + 1] - loglik.seq[1:p]), npara.seq[p + 1] - npara.seq[1:p], lower.tail = F)
  
  if (any(lrt.test > alpha)) {
    u.lrt <- which(lrt.test > alpha)[1] - 1
  } else {
    u.lrt <- p
  }
  
  return(list(u.aic = u.aic, u.bic = u.bic, u.lrt = u.lrt, loglik.seq = loglik.seq, aic.seq = aic.seq, bic.seq = bic.seq))
}

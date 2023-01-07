u.env.apweights <- function(X, Y, alpha = 0.01) {
  
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  loglik.seq <- unlist(lapply(0:r, function(x) env.apweights(X, Y, x, asy = FALSE)$loglik))
  npara.seq <- r + r * (r + 1) / 2 + p * (0:r)
  
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

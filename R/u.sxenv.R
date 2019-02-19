u.sxenv <- function(X, Y, R){
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  q <- length(R)
  loglik.seq <- unlist(lapply(0:p, function(x) sxenv(X, 
                                                Y, x, R, asy = F)$loglik))
  npara.seq <- r * (0:p) + r + p + q - 1 + p * (p + 1)/2 + r * (r + 1)/2
  aic.seq <- -2 * loglik.seq + 2 * npara.seq
  bic.seq <- -2 * loglik.seq + log(n) * npara.seq
  u.aic <- which.min(aic.seq) - 1
  u.bic <- which.min(bic.seq) - 1
  if (min(aic.seq) == aic.seq[p+1]) {
      u.aic = p
      u.bic = p
  }
  return(list(u.aic = u.aic, u.bic = u.bic, aic.seq = aic.seq, 
              bic.seq = bic.seq))
}
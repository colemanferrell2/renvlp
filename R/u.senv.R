u.senv <- function(X, Y){
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  loglik.seq <- unlist(lapply(0:r, function(x) senv(X, Y, x, 
                                                   asy = F)$loglik))
  npara.seq <- 2 * r - 1 + r * (r + 1)/2 + p * (0:r)
  aic.seq <- -2 * loglik.seq + 2 * npara.seq
  bic.seq <- -2 * loglik.seq + log(n) * npara.seq
  u.aic <- which.min(aic.seq) - 1
  u.bic <- which.min(bic.seq) - 1
  if (min(aic.seq) == aic.seq[r+1]) {
      u.aic = r
      u.bic = r
  }
  return(list(u.aic = u.aic, u.bic = u.bic, 
              aic.seq = aic.seq, bic.seq = bic.seq))
}
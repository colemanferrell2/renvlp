u.eppls <- function(X1, X2, Y) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  loglik.seq <- unlist(lapply(0:p1, function(x) eppls(X1, X2, Y, x, asy = FALSE)$loglik))
  npara.seq <- r * (0:p1) + (0:p1) * (1:(p1+1)) / 2 + (p1:0) * ((p1+1):1) / 2 + (0:p1)*(p1:0) + r*(r+1)/2 + r + p1 + p2 + p2 * r + p1 * p2
  bic.seq <- -2 * loglik.seq + log(n) * npara.seq
  u.bic <- which.min(bic.seq) - 1
  return(list(u.bic = u.bic, bic.seq = bic.seq))
}


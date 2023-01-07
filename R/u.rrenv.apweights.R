u.rrenv.apweights <- function(X, Y, d, alpha = 0.01) {
  
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  likl <- rep(0, r - d + 1)
  aic <- rep(0, r - d + 1)
  bic <- rep(0, r - d + 1)
  Npara <- rep(0, r - d + 1)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)

  for (i in d:r) {
    tmp <- rrenv.apweights(X, Y, i, d, asy = FALSE)
    likl[i+1-d] <- tmp$loglik	
    Npara[i+1-d] <- (p + i - d) * d + r * (r + 1) / 2		
    aic[i+1-d] <- 2 * Npara[i+1-d] - 2 * likl[i+1-d]
    bic[i+1-d] <- log(n) * Npara[i+1-d] - 2 * likl[i+1-d]
  }
  
  uhat.aic <- which.min(aic) + d - 1
  uhat.bic <- which.min(bic) + d - 1
  
  lrt.test <- stats::pchisq(2 * (likl[r - d + 1] - likl[1:(r-d)]), Npara[r - d + 1] - Npara[1:(r-d)], lower.tail = F)
  
  if (any(lrt.test > alpha)) {
    u.lrt <- which(lrt.test > alpha)[1] + d - 1
  } else {
    u.lrt <- r
  }
  
  return(list(u.aic = uhat.aic, u.bic = uhat.bic, u.lrt = u.lrt, loglik.seq = likl, aic.seq = aic, bic.seq = bic))
}




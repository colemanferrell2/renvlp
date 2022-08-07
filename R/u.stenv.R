u.stenv <- function(X, Y, alpha = 0.01) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  tmp <- d.stenv(X, Y)
  d <- tmp$rank.beta
  if (d == 0) 
    stop("d must be an integer between 0 and min(r, p).")
  loglik.mat <- matrix(rep(0, p * r), nrow = p)
  for (i in 1 : p) {
    loglik.mat[i , ] <- unlist(lapply(1 : r, function(x) stenv(X, Y, i, x,  
                                                               asy = F)$loglik))
  }
  l.full <- loglik.mat[p, r]
  loglik.seqdx <- loglik.mat[ , d]
  npara.seqdx <- p * r - (d : (p - 1)) * d 
  lrt.testdx <- stats::pchisq(2 * (l.full -loglik.seqdx[d : (p - 1)]),
                       npara.seqdx[1 : (p - d)], lower.tail = F)
  if (any(lrt.testdx > alpha)) {
    u.lrtdx <- which(lrt.testdx > alpha)[1] - 1 + d
  }
  else {
    u.lrtdx <- p
  }
  q <- u.lrtdx
  loglik.seqdy <- loglik.mat[q , ]
  npara.seqdy <- p * r - q * (d : (r - 1)) 
  lrt.testdy <- stats::pchisq(2 * (l.full -loglik.seqdy[d : (r - 1)]), 
                       npara.seqdy[1 : (r - d)], lower.tail = F)
  if (any(lrt.testdy > alpha)) {
    u.lrtdy <- which(lrt.testdy > alpha)[1] - 1 + d
  }
  else {
    u.lrtdy <- r
  }
  u <- u.lrtdy
  u.lrt <- c(q, u)
  npara.mat <- matrix(rep(0, p * r), nrow = p)
  for(i in 1 : p){
    for(j in 1 : r){
      npara.mat[i, j] <- p * (p + 3) / 2 + r * (r + 3) / 2 + i * j
    }
  }
  aic.mat <- -2 * loglik.mat + 2 * npara.mat
  bic.mat <- -2 * loglik.mat + log(n) * npara.mat
  u.aic <- arrayInd(which.min(aic.mat[d : p, d : r]), 
                    dim(as.matrix(aic.mat[d : p, d : r]))) + d - 1
  u.bic <- arrayInd(which.min(bic.mat[d : p, d : r]), 
                    dim(as.matrix(bic.mat[d : p, d : r]))) + d - 1
  
  return(list(d = d, u.aic = u.aic, u.bic = u.bic, u.lrt = u.lrt, 
              loglik.mat = loglik.mat, aic.mat = aic.mat, bic.mat = bic.mat))
}

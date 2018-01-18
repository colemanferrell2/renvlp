u.genv <- function(X, Y, Z, alpha = 0.01) {
  
  ZZ <- as.factor(Z)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(X)
  p <- nlevels(ZZ)
  
  if (n <= c){
    print("This works when p < n")
    skip <- 1
  }else{
    skip <- 0
  }
  
  if(skip == 0){
    
    paranum = loglik.seq <- c()
    for (u in 0 : r){
      loglik.seq[u + 1] <- genv(X, Y, Z, u, asy = F, fit = F)$loglik
      paranum[u + 1] <- p * r + p * u * c + u * (r - u) + p * u * (u + 1)/2 + (r - u)*(r - u + 1)/2
    }
    aic.seq <- -2 * loglik.seq + 2 * paranum
    bic.seq <- -2 * loglik.seq + log(n) * paranum
    u.aic <- which.min(aic.seq) - 1
    u.bic <- which.min(bic.seq) - 1
    lrt.test <- stats::pchisq(2 * (loglik.seq[r + 1] - loglik.seq[1:r]), 
                       paranum[r + 1] - paranum[1:r], lower.tail = F)
    if (any(lrt.test > alpha)) {
      u.lrt <- which(lrt.test > alpha)[1] - 1
    }
    else {
      u.lrt <- r
    }
    
    return(list(u.aic = u.aic, u.bic = u.bic, u.lrt = u.lrt, 
                loglik.seq = loglik.seq, aic.seq = aic.seq, bic.seq = bic.seq))
  }
}
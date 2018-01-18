boot.henv <- function(X, Y, u, B) {

  XX <- as.factor(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- nlevels(XX)
  if(u < 0 | u > r){
    stop("u should be an interger between 0 and r")
  } 
  
  ncumx <- c()
  for (i in 1:p) {
    ncumx[i] <- length(which(XX == as.numeric(levels(XX)[i])))
  }
  ncum <- cumsum(ncumx)
  ng <- diff(c(0, ncum))
  sortx <- sort(X, index.return = T)
  Xs <- sortx$x
  ind <- sortx$ix
  fit <- henv(X, Y, u, asy = F, fit = T)
  Yfit <- fit$Yfit
  res <- Y - Yfit
  boothenv <- function(i) {
    res.boot <- matrix(rep(0, n * r), ncol = r)
    for (j in 1:p) {
      if (j > 1) {
        res.boot[ind[(ncum[j - 1] + 1):ncum[j]], ] <- as.matrix(res[sample(ind[(ncum[j - 1] + 1)
                                                      :ncum[j]], ng[j], replace = T), ], ncol = r)
      }
      else {
        res.boot[ind[1:ncum[1]], ] <- as.matrix(res[sample(ind[1:ncum[1]], ng[1], replace = T), ], ncol = r)
      }
    }
    Y.boot <- Yfit + res.boot
    return(c(henv(X, Y.boot, u, asy = F)$beta))
  }
  bootbeta <- lapply(1:B, function(i) boothenv(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)
  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = r)
  return(bootse)
}
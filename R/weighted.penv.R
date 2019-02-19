weighted.penv <- function(X1, X2, Y, bstrpNum = 0, min.u = 1, max.u = ncol(as.matrix(Y)), boot.resi = "full") {
  
  bstrpNum <- ceiling(bstrpNum)
  if (bstrpNum < 0) stop("The number of bootstrap samples should be a positive integer.")
  if (min.u < 0) stop("The smallest dimension of the envelope subspace is 0.")
  if (max.u > ncol(Y)) stop("The largest dimension of the envelope subspace is r.")
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Y <- as.matrix(Y)
  bootse <- NULL
  ratios <- NULL
  bic_select <- NULL
  
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p1 <- ncol(X1) 
  p2 <- ncol(X2) 
  p <- p1 + p2
  
  w.penv <- function(X1, X2, Y, min.u, max.u) {
  
    min.u <- ceiling(min.u)
    max.u <- ceiling(max.u)
    beta1array <- array(0, c(r, p1, max.u-min.u+1))
    beta2array <- array(0, c(r, p2, max.u-min.u+1))
    muarray <- matrix(0, r, max.u-min.u+1)
    Sigmaarray <- array(0, c(r, r, max.u-min.u+1))
    bic.seq <- rep(0, max.u-min.u+1)
    w <- rep(0, max.u-min.u+1)
    beta1 <- matrix(0, r, p1)
    beta2 <- matrix(0, r, p2)
    mu <- matrix(0, r, 1)
    Sigma <- matrix(0, r, r)
    
    for (i in min.u:max.u) {
      m <- penv(X1, X2, Y, i, asy = F)
      beta1array[ , , i-min.u+1] <- m$beta1
      beta2array[ , , i-min.u+1] <- m$beta2
      muarray[, i-min.u+1] <- m$mu
      Sigmaarray[ , , i-min.u+1] <- m$Sigma
      bic.seq[i-min.u+1] <- -2 * m$loglik + log(n) * (r + r * p2 + r * (r + 1) / 2 + p1 * i)
    }
    
    for (i in min.u:max.u) {
      w[i-min.u+1] <- 1 / sum(exp(bic.seq[i-min.u+1] - bic.seq))
      beta1 <- w[i-min.u+1] * beta1array[ , , i-min.u+1] + beta1
      beta2 <- w[i-min.u+1] * beta2array[ , , i-min.u+1] + beta2
      mu <- w[i-min.u+1] * muarray[ , i-min.u+1] + mu
      Sigma <- w[i-min.u+1] * Sigmaarray[ , , i-min.u+1] + Sigma
    }
    
    return(list(beta1 = beta1, beta2 = beta2, mu = mu, w = w, Sigma = Sigma))
  }
  
  m <- w.penv(X1, X2, Y, min.u, max.u)
  mu <- m$mu  
  beta1 <- m$beta1
  beta2 <- m$beta2
  Sigma <- m$Sigma
  w <- m$w
  Yfit <- matrix(1, n, 1) %*% t(mu) + X1 %*% t(beta1) + X2 %*% t(beta2)
  if (boot.resi == "weighted") {
    res <- Y - Yfit
  } else {
    m2 <- penv(X1, X2, Y, r, asy = F)
    mu.full <- m2$mu  
    beta1.full <- m2$beta1
    beta2.full <- m2$beta2
    res <- Y - matrix(1, n, 1) %*% t(mu.full) - X1 %*% t(beta1.full) - X2 %*% t(beta2.full)
  }
  res1 <- as.matrix(res)
  tmp <- res1 %*% chol2inv(chol(Sigma))
  loglik <- - n * r / 2 * log(2 * pi) - n / 2 * sum(log(eigen(Sigma)$values)) - 1 / 2 * sum(c(tmp) * c(res1))
  
  if (bstrpNum > 0) {
    bootpenv <- function(i) {
      res.boot <- res[sample(1:n, n, replace = T), ]
      Y.boot <- Yfit + res.boot
      return(c(w.penv(X1, X2, Y.boot, min.u, max.u)$beta1, w.penv(X1, X2, Y.boot, min.u, max.u)$w))
    }
    
    bootres <- lapply(1:bstrpNum, function(i) bootpenv(i))
    bootres <- matrix(unlist(bootres), nrow = bstrpNum, byrow = TRUE)
    
    bootse <- matrix(apply(bootres[, 1:(r*p1)], 2, stats::sd), nrow = r)
    bic_select <- table(apply(as.matrix(bootres[, (r*p1+1):(r*p1+max.u-min.u+1)]), 1, which.max))
    names(bic_select) <- as.character(strtoi(names(bic_select)) + min.u - 1)
    
    bootse.full <- boot.penv(X1, X2, Y, r, bstrpNum)
    ratios <- bootse.full / bootse
    
  }
  
  return(list(beta1 = beta1, beta2 = beta2, mu = mu, w = w, Sigma = Sigma, loglik = loglik, n = n, bootse = bootse, ratios = ratios, bic_select = bic_select))
}

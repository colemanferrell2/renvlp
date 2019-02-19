weighted.xenv <- function(X, Y, bstrpNum = 0, min.u = 1, max.u = ncol(as.matrix(X)), boot.resi = "full") {
  
  bstrpNum <- ceiling(bstrpNum)
  if (bstrpNum < 0) stop("The number of bootstrap samples should be a positive integer.")
  if (min.u < 0) stop("The smallest dimension of the envelope subspace is 0.")
  if (max.u > ncol(X)) stop("The largest dimension of the envelope subspace is p.")
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  bootse <- NULL
  ratios <- NULL
  bic_select <- NULL
  
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X) 
  
  w.xenv <- function(X, Y, min.u, max.u) {
    
    min.u <- ceiling(min.u)
    max.u <- ceiling(max.u)
    betaarray <- array(0, c(p, r, max.u-min.u+1))
    muarray <- matrix(0, r, max.u-min.u+1)
    SigmaXarray <- array(0, c(p, p, max.u-min.u+1))
    SigmaYcXarray <- array(0, c(r, r, max.u-min.u+1))
    bic.seq <- rep(0, max.u-min.u+1)
    w <- rep(0, max.u-min.u+1)
    beta <- matrix(0, p, r)
    mu <- matrix(0, r, 1)
    SigmaX <- matrix(0, p, p)
    SigmaYcX <- matrix(0, r, r)
      
    for (i in min.u:max.u) {
      m <- xenv(X, Y, i, asy = F)
      betaarray[ , , i-min.u+1] <- m$beta
      muarray[, i-min.u+1] <- m$mu
      SigmaXarray[ , , i-min.u+1] <- m$SigmaX
      SigmaYcXarray[ , , i-min.u+1] <- m$SigmaYcX
      bic.seq[i-min.u+1] <- -2 * m$loglik + log(n) * (r + r * (r + 1) / 2 + p * (p + 1) / 2 + r * i)
    }
    
    for (i in min.u:max.u) {
      w[i-min.u+1] <- 1 / sum(exp(bic.seq[i-min.u+1] - bic.seq))
      beta <- w[i-min.u+1] * betaarray[ , , i-min.u+1] + beta
      mu <- w[i-min.u+1] * muarray[ , i-min.u+1] + mu
      SigmaX <- w[i-min.u+1] * SigmaXarray[ , , i-min.u+1] + SigmaX
      SigmaYcX <- w[i-min.u+1] * SigmaYcXarray[ , , i-min.u+1] + SigmaYcX
    }
    
    return(list(beta = beta, mu = mu, w = w, SigmaX = SigmaX, SigmaYcX = SigmaYcX))
  }
  
  m <- w.xenv(X, Y, min.u, max.u)
  beta <- m$beta
  mu <- m$mu
  SigmaX <- m$SigmaX
  SigmaYcX <- m$SigmaYcX
  w <- m$w
  Yfit <- matrix(1, n, 1) %*% t(mu) + X %*% beta
  if (boot.resi == "weighted") {
    res <- Y - Yfit
  } else {
    m2 <- xenv(X, Y, p, asy = F)
    beta.full <- m2$beta
    mu.full <- m2$mu
    res <- Y - matrix(1, n, 1) %*% t(mu.full) - X %*% beta.full
  }
  res1 <- as.matrix(res)
  tmp <- res1 %*% chol2inv(chol(SigmaYcX))
  resX <- scale(X, center = T, scale = F)
  tmpX <- resX %*% chol2inv(chol(SigmaX))
  loglik <- - n * (r + p) / 2 * log(2 * pi) - n / 2 * sum(log(eigen(SigmaYcX)$values)) - 1 / 2 * sum(c(tmp) * c(res1)) - n / 2 * sum(log(eigen(SigmaX)$values)) - 1 / 2 * sum(c(tmpX) * c(resX))
    
  if (bstrpNum > 0) {
    bootenv <- function(i) {
      res.boot <- res[sample(1:n, n, replace = T), ]
      Y.boot <- Yfit + res.boot
      return(c(w.xenv(X, Y.boot, min.u, max.u)$beta, w.xenv(X, Y.boot, min.u, max.u)$w))
    }
    
    bootres <- lapply(1:bstrpNum, function(i) bootenv(i))
    bootres <- matrix(unlist(bootres), nrow = bstrpNum, byrow = TRUE)
    
    bootse <- matrix(apply(bootres[, 1:(p*r)], 2, stats::sd), nrow = p)
    bic_select <- table(apply(as.matrix(bootres[, (r*p+1):(r*p+max.u-min.u+1)]), 1, which.max))
    names(bic_select) <- as.character(strtoi(names(bic_select)) + min.u - 1)
    
    bootse.full <- boot.xenv(X, Y, p, bstrpNum)
    ratios <- bootse.full / bootse
    
  }
  
  return(list(beta = beta, mu = mu, SigmaX = SigmaX, SigmaYcX = SigmaYcX, w = w, loglik = loglik, n = n, bootse = bootse, ratios = ratios, bic_select = bic_select))
}

weighted.env <- function(X, Y, bstrpNum = 0, min.u = 1, max.u = ncol(Y), boot.resi = "full") {
  
  bstrpNum <- ceiling(bstrpNum)
  if (bstrpNum < 0) stop("The number of bootstrap samples should be a positive integer.")
  if (min.u < 0) stop("The smallest dimension of the envelope subspace is 0.")
  if (max.u > ncol(Y)) stop("The largest dimension of the envelope subspace is r.")
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  bootse <- NULL
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)  
  select.seq <- NULL
  ratios <- NULL
  bic_select <- NULL
  
  w.env <- function(X, Y, min.u, max.u) {
 
    min.u <- ceiling(min.u)
    max.u <- ceiling(max.u)
    betaarray <- array(0, c(r, p, max.u-min.u+1))
    muarray <- matrix(0, r, max.u-min.u+1)
    Sigmaarray <- array(0, c(r, r, max.u-min.u+1))
    bic.seq <- rep(0, max.u-min.u+1)
    w <- rep(0, max.u-min.u+1)
    beta <- matrix(0, r, p)
    mu <- matrix(0, r, 1)
    Sigma <- matrix(0, r, r)
      
    for (i in min.u:max.u) {
      m <- env(X, Y, i, asy = F)
      betaarray[ , , i-min.u+1] <- m$beta
      muarray[, i-min.u+1] <- m$mu
      Sigmaarray[ , , i-min.u+1] <- m$Sigma
      bic.seq[i-min.u+1] <- -2 * m$loglik + log(n) * (r + r * (r + 1) / 2 + p * i)
    }
    
    for (i in min.u:max.u) {
      w[i-min.u+1] <- 1 / sum(exp(bic.seq[i-min.u+1] - bic.seq))
      beta <- w[i-min.u+1] * betaarray[ , , i-min.u+1] + beta
      mu <- w[i-min.u+1] * muarray[ , i-min.u+1] + mu
      Sigma <- w[i-min.u+1] * Sigmaarray[ , , i-min.u+1] + Sigma
    }
    
    return(list(beta = beta, mu = mu, w = w, Sigma = Sigma))
  }
  
  m <- w.env(X, Y, min.u, max.u)
  mu <- m$mu  
  beta <- m$beta
  Sigma <- m$Sigma
  w <- m$w
  Yfit <- matrix(1, n, 1) %*% t(mu) + X %*% t(beta)
  if (boot.resi == "weighted") {
    res <- Y - Yfit
  } else {
    m2 <- env(X, Y, r, asy = F)
    mu.full <- m2$mu  
    beta.full <- m2$beta
    res <- Y - matrix(1, n, 1) %*% t(mu.full) - X %*% t(beta.full)
  }
  res1 <- as.matrix(res)
  tmp <- res1 %*% chol2inv(chol(Sigma))
  loglik <- - n * r / 2 * log(2 * pi) - n / 2 * sum(log(eigen(Sigma)$values)) - 1 / 2 * sum(c(tmp) * c(res1))
    
  if (bstrpNum > 0) {
    bootenv <- function(i) {
      res.boot <- res[sample(1:n, n, replace = T), ]
      Y.boot <- Yfit + res.boot
      return(c(w.env(X, Y.boot, min.u, max.u)$beta, w.env(X, Y.boot, min.u, max.u)$w))
    }
    
    bootres <- lapply(1:bstrpNum, function(i) bootenv(i))
    bootres <- matrix(unlist(bootres), nrow = bstrpNum, byrow = TRUE)
    
    bootse <- matrix(apply(bootres[, 1:(r*p)], 2, stats::sd), nrow = r)
    bic_select <- table(apply(as.matrix(bootres[, (r*p+1):(r*p+max.u-min.u+1)]), 1, which.max))
    names(bic_select) <- as.character(strtoi(names(bic_select)) + min.u - 1)
    
    
    bootse.full <- boot.env(X, Y, r, bstrpNum)
    ratios <- bootse.full / bootse
  }
  
  return(list(beta = beta, mu = mu, Sigma = Sigma, w = w, loglik = loglik, n = n, bootse = bootse, ratios = ratios, bic_select = bic_select))
}

penv <- function(X1, X2, Y, u, asy = TRUE, init = NULL) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  
  if (a[1] != nrow(X1)) stop("X1 and Y should have the same number of observations.")
  if (a[1] != nrow(X2)) stop("X2 and Y should have the same number of observations.")
  if (u > r || u < 0) stop("u must be an integer between 0 and r.")
  if(sum(duplicated(cbind(X1, X2, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
      if (nrow(init) != r || ncol(init) != u) stop("The dimension of init is wrong.")
  }

  Yc <- scale(Y, center = T, scale = F)
  X1c <- scale(X1, center = T, scale = F)
  X2c <- scale(X2, center = T, scale = F)
  sigY <- stats::cov(Y) * (n - 1) / n
  sigX2 <- stats::cov(X2) * (n - 1) / n
  invsigX2 <- chol2inv(chol(sigX2))
  sigX2X1 <- stats::cov(X2, X1) * (n - 1) / n
  sigX2Y <- stats::cov(X2, Y) * (n - 1) / n

  res.1c2 <- X1c - X2c %*% invsigX2 %*% sigX2X1
  res.yc2 <- Yc - X2c %*% invsigX2 %*% sigX2Y


  tmp <- env(res.1c2, res.yc2, u, init = init)
  
  
  Gammahat <- tmp$Gamma
  Gamma0hat <- tmp$Gamma0
  loglik <- tmp$loglik
  
  covMatrix <- NULL
  asySE1 <- NULL
  asySE2 <- NULL
  ratio <- NULL
  
  if (u == 0) {
    
    beta1hat <- matrix(0, r, p1)
    beta2hat <- crossprod(sigX2Y, invsigX2)
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- tmp$Sigma
    muhat <- colMeans(Y) - beta2hat %*% colMeans(X2)
    Sigmahat <- Omega0hat
    if (asy == T) {
      covMatrix <- kronecker(invsigX2, Sigmahat)
      asySE1 <- matrix(0, r, p1)
      asySE2 <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- matrix(1, r, p1)
    }
    
  } else if (u == r) {
    
    p <- p1 + p2
    X <- cbind(X1, X2)
    sigYX <- stats::cov(Y, X) * (n - 1) / n
    sigX <- stats::cov(X) * (n - 1) / n
    invsigX <- chol2inv(chol(sigX))
    betaOLS <- sigYX %*% invsigX
    beta1hat <- betaOLS[ , 1:p1]
    beta2hat <- betaOLS[ , (p1 + 1):p]
    etahat <- beta1hat
    Omegahat <- tmp$Sigma
    Omega0hat <- NULL
    muhat <- colMeans(Y) - betaOLS %*% colMeans(X)
    Sigmahat <- Omegahat
    if (asy == T) {
      covMatrix <- kronecker(invsigX, Sigmahat)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      asySE1 <- asySE[, 1:p1, drop = F]
      asySE2 <- asySE[, (p1+1):p]
      ratio <- matrix(1, r, p1)
    }
    
  } else {
    
    etahat <- tmp$eta
    beta1hat <- tmp$beta
    beta2hat <- (t(sigX2Y) - tcrossprod(beta1hat, sigX2X1)) %*% invsigX2
    muhat <- colMeans(Y) - beta1hat %*% colMeans(X1) - beta2hat %*% colMeans(X2)
    Omegahat <- tmp$Omega
    Omega0hat <- tmp$Omega0
    Sigmahat <- tmp$Sigma
    if (asy == T) {
      sig.1c2 <- stats::cov(res.1c2) * (n - 1) / n
      invsig.1c2 <- chol2inv(chol(sig.1c2))
      sig.yc2 <- stats::cov(res.yc2) * (n - 1) / n
      sig.yc2.1c2 <- stats::cov(res.yc2, res.1c2) * (n - 1) / n
      sig.ycx <- sig.yc2 - sig.yc2.1c2 %*% tcrossprod(invsig.1c2, sig.yc2.1c2)
      covMatrix <- kronecker(invsig.1c2, sig.ycx)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
      
      invOmega0hat <- chol2inv(chol(Omega0hat))
      invOmegahat <- chol2inv(chol(Omegahat))
      temp <- kronecker(etahat %*% tcrossprod(sig.1c2, etahat), invOmega0hat) + kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u))
      temp2 <- kronecker(t(etahat), Gamma0hat)
      Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
      covMatrix <- kronecker(invsig.1c2, Sigma1) + temp2 %*% chol2inv(chol(temp)) %*% t(temp2)
      asySE1 <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- asyFm / asySE1
      
      p <- p1 + p2
      invsig.ycx <- chol2inv(chol(sig.ycx))
      sigX1 <- stats::cov(X1) * (n - 1) / n
      sep1 <- p2 * r
      sep2 <- p * r
      sep3 <- r * p2 + u * p1
      sep4 <- r * p2 + u * p1 + u * (r - u)
      sep5 <- r * p2 + u * p1 + u * (r - u) + u * (u + 1) / 2
      J <- diag(r * p + r * (r + 1) / 2)
      J[1:sep1, 1:sep1] <- kronecker(sigX2, invsig.ycx)
      J[1:sep1, (sep1 + 1):sep2] <- kronecker(sigX2X1, invsig.ycx)
      J[(sep1 + 1):sep2, 1:sep1] <- t(J[1:sep1, (sep1 + 1):sep2])
      J[(sep1 + 1):sep2, (sep1 + 1):sep2] <- kronecker(sigX1, invsig.ycx)
      J[(sep2 + 1):nrow(J), (sep2 + 1):ncol(J)] <- 0.5 * crossprod(expan(r), kronecker(invsig.ycx, invsig.ycx)) %*% expan(r) 
      H <- matrix(0, r * p + r * (r + 1) / 2, r * p2 + u * p1 + r * (r + 1) / 2)
      H[1:sep1, 1:sep1] <- diag(p2 * r)
      H[(sep1 + 1):sep2, (sep1 + 1):sep3] <- kronecker(diag(p1), Gammahat)
      H[(sep1 + 1):sep2, (sep3 + 1):sep4] <- temp2
      H[(sep2 + 1):nrow(H), (sep3 + 1):sep4] <- 2 * contr(r) %*% (kronecker(Gammahat %*% Omegahat, Gamma0hat) - kronecker(Gammahat, Gamma0hat %*% Omega0hat))
      H[(sep2 + 1):nrow(H), (sep4 + 1):sep5] <- contr(r) %*% kronecker(Gammahat, Gammahat) %*% expan(u)
      H[(sep2 + 1):nrow(H), (sep5 + 1):ncol(H)] <- contr(r) %*% kronecker(Gamma0hat, Gamma0hat) %*% expan(r - u)
      
      covMatrix <- H %*% tcrossprod(chol2inv(chol(crossprod(H, J) %*% H)), H)
      covMatrix <- covMatrix[1:sep2, 1:sep2]
      covMatrix <- covMatrix[c((sep1 + 1):nrow(covMatrix), 1:sep1), c((sep1 + 1):nrow(covMatrix), 1:sep1)]
      asySE2 <- matrix(sqrt(diag(covMatrix[(p1*r+1):(p*r), (p1*r+1):(p*r)])), nrow = r)
    }
  }
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat, beta1 = beta1hat, beta2 = beta2hat, Sigma = Sigmahat, eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik, n = n, covMatrix = covMatrix, asySE1 = asySE1, asySE2 = asySE2, ratio = ratio))
}


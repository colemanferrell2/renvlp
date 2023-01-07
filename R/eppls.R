eppls <- function(X1, X2, Y, u, asy = TRUE, init = NULL) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  p <- p1 + p2
  
  if (a[1] != nrow(X1)) stop("X1 and Y should have the same number of observations.")
  if (a[1] != nrow(X2)) stop("X2 and Y should have the same number of observations.")
  if (u > p1 || u < 0) stop("u must be an integer between 0 and p1.")
  if (sum(duplicated(cbind(X1, X2, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
      if (nrow(init) != p1 || ncol(init) != u) stop("The dimension of init is wrong.")
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
  sig.res.1c2 <- stats::cov(res.1c2) * (n - 1) / n
  sig.res.yc2 <- stats::cov(res.yc2) * (n - 1) / n
  sig.res.yx1 <- stats::cov(res.yc2, res.1c2) * (n - 1) / n
  
  if (!is.null(init)) {
    tmp <- xenv(res.1c2, res.yc2, u, asy = FALSE, init = init)
  } else if (n < 100 & u != 0) {
    Xc <- scale(res.1c2, center = T, scale = F)
    Ycr <- scale(res.yc2, center = T, scale = F)
    respls <- pls::mvr(Ycr ~ Xc, center = TRUE, ncomp = u, validation = "CV", method = "simpls")
    initmat <- respls$loadings
    tmp <- xenv(res.1c2, res.yc2, u, init = initmat)
  } else {
    tmp <- xenv(res.1c2, res.yc2, u)
  }
  
  Gammahat <- tmp$Gamma
  Gamma0hat <- tmp$Gamma0
  
  covMatrix1 <- tmp$covMatrix
  
  covMatrix2 <- NULL
  asySE1 <- NULL
  asySE2 <- NULL
  
  if (u == 0) {
    
    gammahat <- crossprod(invsigX2, sigX2X1)
    beta1hat <- matrix(0, p1, r)
    beta2hat <- crossprod(invsigX2, sigX2Y)
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- tmp$SigmaX
    muYhat <- colMeans(Y)
    mu1hat <- colMeans(X1)
    mu2hat <- colMeans(X2)
    SigmaX1hat <- Omega0hat
    SigmaYcXhat <- sig.res.yc2 
    
    eigtemOmega0hat <- eigen(Omega0hat)$values
    logDetOmega0hat <- log(prod(eigtemOmega0hat[eigtemOmega0hat > 0]))
    eigtemSigmaYcXhat <- eigen(SigmaYcXhat)$values
    logDetSigYXhat <- log(prod(eigtemSigmaYcXhat[eigtemSigmaYcXhat > 0]))
    traux <- (X1c - X2c %*% gammahat) %*% solve(SigmaX1hat) %*% t(X1c - X2c %*% gammahat)+
      (Yc - X1c %*% beta1hat - X2c %*% beta2hat) %*% solve(SigmaYcXhat) %*% t(Yc - X1c %*% beta1hat - X2c %*% beta2hat)
    loglik <- -n * (r + p1)/2 * log(2 * pi)-n * logDetOmega0hat / 2 - n * logDetSigYXhat / 2 - sum(diag(traux)) / 2
    if (asy == T) {
        covMatrix2 <- kronecker(SigmaYcXhat, invsigX2)
        asySE2 <- matrix(sqrt(diag(covMatrix2)), nrow = p2)
    }
    
  } else if (u == p1) {
    
    gammahat <- crossprod(invsigX2, sigX2X1)
    etaaux <- crossprod(Gammahat, sig.res.1c2) %*% Gammahat
    invetaaux <- chol2inv(chol(etaaux))
    etaaux2 <- tcrossprod(Gammahat, sig.res.yx1)
    etahat <-  invetaaux %*% etaaux2
    Omegahat <- tmp$SigmaX
    Omega0hat <- NULL
    beta1hat <- tmp$beta
    beta2hat <- invsigX2 %*% (sigX2Y - sigX2X1 %*% beta1hat)
    muYhat <- colMeans(Y)
    mu1hat <- colMeans(X1)
    mu2hat <- colMeans(X2)
    SigmaX1hat <- Omegahat
    SigmaYcXaux <- sig.res.yx1 %*% Gammahat
    SigmaYcXhat <- sig.res.yc2 - SigmaYcXaux %*% invetaaux %*% t(SigmaYcXaux)
    
    eigtemOmegahat <- eigen(Omegahat)$values
    logDetOmegahat <- log(prod(eigtemOmegahat[eigtemOmegahat > 0]))
    eigtemSigmaYcXhat <- eigen(SigmaYcXhat)$values
    logDetSigYXhat <- log(prod(eigtemSigmaYcXhat[eigtemSigmaYcXhat > 0]))
    
    if (matrixcalc::is.singular.matrix(SigmaYcXhat, tol = 1e-08) == F) {
        ssy <- solve(SigmaYcXhat)
    } else {
        ssy <- ginv(SigmaYcXhat, tol = 1e-08)
    }
    
    traux <- (X1c - X2c %*% gammahat) %*% solve(SigmaX1hat) %*% t(X1c - X2c %*% gammahat) +
        (Yc - X1c %*% beta1hat - X2c %*% beta2hat) %*% ssy %*% t(Yc - X1c %*% beta1hat - X2c %*% beta2hat)
    loglik <- -n * (r + p1)/2 * log(2 * pi) - n * logDetOmegahat / 2 - n * logDetSigYXhat / 2 - sum(diag(traux)) / 2
        
    if (asy == T) {
        invOmegahat <- chol2inv(chol(Omegahat))
        covMatrix2 <- kronecker(SigmaYcXhat,invsigX2) + kronecker(SigmaYcXhat, gammahat %*% invOmegahat %*% t(gammahat))
    
        asySE2 <- matrix(sqrt(diag(covMatrix2)), nrow = p2)
    }
    
  } else {
    
    gammahat <- crossprod(invsigX2, sigX2X1)
    etaaux <- crossprod(Gammahat, sig.res.1c2) %*% Gammahat
    invetaaux <- chol2inv(chol(etaaux))
    etaaux2 <- crossprod(Gammahat, t(sig.res.yx1))
    etahat <- invetaaux %*% etaaux2
    Omegahat <- tmp$Omega
    Omega0hat <- tmp$Omega0
    beta1hat <- tmp$beta
    beta2hat <- invsigX2 %*% (sigX2Y - sigX2X1 %*% beta1hat)
    muYhat <- colMeans(Y)
    mu1hat <- colMeans(X1)
    mu2hat <- colMeans(X2)
    SigmaX1hat <- tmp$SigmaX 
    
    SigmaYcXaux <- sig.res.yx1 %*% Gammahat
    SigmaYcXhat <- sig.res.yc2 - SigmaYcXaux %*% invetaaux %*% t(SigmaYcXaux)
    
    eigtemOmegahat <- eigen(Omegahat)$values
    logDetOmegahat <- log(prod(eigtemOmegahat[eigtemOmegahat > 0]))
    eigtemOmega0hat <- eigen(Omega0hat)$values
    logDetOmega0hat <- log(prod(eigtemOmega0hat[eigtemOmega0hat > 0]))
    eigtemSigmaYcXhat <- eigen(SigmaYcXhat)$values
    logDetSigYXhat <- log(prod(eigtemSigmaYcXhat[eigtemSigmaYcXhat > 0]))
    
    if (matrixcalc::is.singular.matrix(SigmaYcXhat, tol = 1e-08) == F) {
        ssy <- solve(SigmaYcXhat)
    } else {
        ssy <- ginv(SigmaYcXhat, tol = 1e-08)
    }
    
    traux <- (X1c - X2c %*% gammahat) %*% solve(SigmaX1hat) %*% t(X1c - X2c %*% gammahat) +
      (Yc - X1c %*% beta1hat - X2c %*% beta2hat) %*% ssy %*% t(Yc - X1c %*% beta1hat - X2c %*% beta2hat)
    loglik <- -n * (r + p1)/2 * log(2 * pi) - n * logDetOmegahat / 2 - n * logDetOmega0hat / 2 - n * logDetSigYXhat / 2 - sum(diag(traux)) / 2
    
    if (asy == T) {
        invOmegahat <- chol2inv(chol(Omegahat))
        invOmega0hat <- chol2inv(chol(Omega0hat))
        tempM <- kronecker(etahat %*% tcrossprod(ssy, etahat), Omega0hat) + kronecker(invOmegahat,
            Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(p1 - u))
        invMmd <- chol2inv(chol(tempM))
        temp2 <- gammahat %*% Gammahat %*% invOmegahat %*% t(Gammahat) %*% t(gammahat)
        covMatrix2 <- kronecker(SigmaYcXhat, invsigX2) + kronecker(SigmaYcXhat, temp2) +
            kronecker(t(etahat), gammahat %*% Gamma0hat) %*% invMmd %*% kronecker(etahat, t(Gamma0hat) %*% t(gammahat))
        asySE2 <- matrix(sqrt(diag(covMatrix2)), nrow = p2)
    }
    
  }
  
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, gamma = gammahat,
              muY = muYhat, mu1 = mu1hat, mu2 = mu2hat,
              beta1 = beta1hat, beta2 = beta2hat, 
              SigmaX1 = SigmaX1hat, SigmaYcX = SigmaYcXhat, eta = etahat, 
              Omega = Omegahat, Omega0 = Omega0hat, 
              loglik = loglik, n = n, covMatrix1 = covMatrix1,
              covMatrix2 = covMatrix2,
              asySE1 = tmp$asySE, asySE2 = asySE2))
}


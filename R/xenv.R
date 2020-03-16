xenv <- function(X, Y, u, asy = TRUE, init = NULL) {
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  if (a[1] != nrow(X)) stop("X and Y should have the same number of observations.")
  if (u > p || u < 0) stop("u must be an integer between 0 and p.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
    if (nrow(init) != p || ncol(init) != u) stop("The dimension of init is wrong.")
  }
  sigY <- stats::cov(Y) * (n - 1) / n
  sigYX <- stats::cov(Y, X) * (n - 1) / n
  sigX <- stats::cov(X) * (n - 1) / n
  invsigY <- chol2inv(chol(sigY))
  eig.sigY <- eigen(sigY)

  
  U <- crossprod(sigYX, invsigY) %*% sigYX 
  M <- sigX - U
  
  tmp <- envMU(M, U, u, initial = init)
  
  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  objfun <- tmp$objfun
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL 
  
  if (u == 0) {
    
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigX
    muhat <- colMeans(Y)
    betahat <- matrix(0, p, r)                                         
    SigmaXhat <- sigX
    SigmaYcXhat <- sigY
    loglik <- - n * (p + r) / 2 * (log(2 * pi) + 1) - n / 2 * (objfun + sum(log(eig.sigY$values)))
    if (asy == T) ratio <- matrix(1, p, r)
    
  } else if (u == p) {
    
    invsigX <- chol2inv(chol(sigX))
    betaOLS <- tcrossprod(invsigX, sigYX)
    etahat <- betaOLS
    Omegahat <- sigX
    Omega0hat <- NULL
    muhat <- colMeans(Y) - crossprod(betaOLS, colMeans(X))
    betahat <- betaOLS
    SigmaXhat <- M + U
    SigmaYcXhat <- sigY - sigYX %*% betaOLS
    loglik <- - n * (r + p) / 2 * (log(2 * pi) + 1) - n / 2 * (objfun + sum(log(eig.sigY$values)))
    if (asy == T) {
      covMatrix <- kronecker(SigmaYcXhat, invsigX)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- matrix(1, p, r)
    }
    
  } else {
    
    invsigX <- chol2inv(chol(sigX))
    etahat <- crossprod(Gammahat, t(sigYX))
    Omegahat <- crossprod(Gammahat, sigX) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigX) %*% Gamma0hat
    invOmegahat <- chol2inv(chol(Omegahat))
    betahat <- Gammahat %*% invOmegahat %*% etahat
    muhat <- colMeans(Y) - crossprod(betahat, colMeans(X))
    SigmaXhat <- Gammahat %*% tcrossprod(Omegahat, Gammahat) + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
    betaOLS <- tcrossprod(invsigX, sigYX)
    PGamma <- tcrossprod(Gammahat)
    SigmaYcXhat <- sigY - sigYX %*% PGamma %*% chol2inv(chol(SigmaXhat)) %*% PGamma %*% t(sigYX)
    loglik <- - n * (r + p) / 2 * (log(2 * pi) + 1) - n / 2 * (objfun + sum(log(eig.sigY$values)))
    if (asy == T) {
      covMatrix <- kronecker(SigmaYcXhat, invsigX)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = p)
      invSigmaYcXhat <- chol2inv(chol(SigmaYcXhat))
      invOmega0hat <- chol2inv(chol(Omega0hat))
      temp <- kronecker(etahat %*% tcrossprod(invSigmaYcXhat, etahat), Omega0hat) + kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(p - u))
      temp2 <- kronecker(t(etahat), Gamma0hat)
      covMatrix <- kronecker(SigmaYcXhat, Gammahat %*% tcrossprod(invOmegahat, Gammahat)) + temp2 %*% chol2inv(chol(temp)) %*% t(temp2)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- asyFm / asySE
    }    
  }
  
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat, beta = as.matrix(betahat), SigmaX = SigmaXhat, eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, SigmaYcX = SigmaYcXhat, loglik = loglik, n = n, covMatrix = covMatrix, asySE = asySE, ratio = ratio))
}




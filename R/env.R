env <- function(X, Y, u, asy = TRUE, init = NULL) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  if (a[1] != nrow(X)) stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) stop("u must be an integer between 0 and r.")
  if(sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
      if (nrow(init) != r || ncol(init) != u) stop("The dimension of init is wrong.")
  }
  
  sigY <- stats::cov(Y) * (n - 1) / n
  sigYX <- stats::cov(Y, X) * (n - 1) / n
  sigX <- stats::cov(X) * (n - 1) / n
  invsigX <- chol2inv(chol(sigX))
  betaOLS <- sigYX %*% invsigX
  
  U <- tcrossprod(betaOLS, sigYX) 
  M <- sigY - U
  
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
		Omega0hat <- sigY
		muhat <- colMeans(Y)
		betahat <- matrix(0, r, p)
		Sigmahat <- sigY
		loglik <- - n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
		if (asy == T) ratio <- matrix(1, r, p)
		
	} else if (u == r) {
	
	  etahat <- betaOLS
	  Omegahat <- M
	  Omega0hat <- NULL
	  muhat <- colMeans(Y) - betaOLS %*% colMeans(X)
	  betahat <- betaOLS
	  Sigmahat <- M
	  loglik <- - n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
	  if (asy == T) {
	    covMatrix <- kronecker(invsigX, M)
	    asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
	    ratio <- matrix(1, r, p)
	  }
	  
  } else {
    
    etahat <- crossprod(Gammahat, betaOLS)
    betahat <- Gammahat %*% etahat
    muhat <- colMeans(Y) - betahat %*% colMeans(X)
    Omegahat <- crossprod(Gammahat, M) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
    Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
    Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
    loglik <- - n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(invsigX, M)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
      invOmega0hat <- chol2inv(chol(Omega0hat))
      invOmegahat <- chol2inv(chol(Omegahat))
      temp <- kronecker(etahat %*% tcrossprod(sigX, etahat), invOmega0hat) + kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u))
      temp2 <- kronecker(t(etahat), Gamma0hat)
      covMatrix <- kronecker(invsigX, Sigma1) + temp2 %*% chol2inv(chol(temp)) %*% t(temp2)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- asyFm / asySE
    }    
  }
 
   return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat, beta = betahat, Sigma = Sigmahat, eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik, n = n, covMatrix = covMatrix, asySE = asySE, ratio = ratio))
}
	



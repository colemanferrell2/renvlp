env.apweights <- function(X, Y, u, asy = TRUE) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  if (a[1] != nrow(X)) stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) stop("u must be an integer between 0 and r.")
  if(sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  
  sigY <- stats::cov(Y) * (n - 1) / n
  sigYX <- stats::cov(Y, X) * (n - 1) / n
  sigX <- stats::cov(X) * (n - 1) / n
  invsigX <- chol2inv(chol(sigX))
  betaOLS <- sigYX %*% invsigX
  
  U <- tcrossprod(betaOLS, sigYX) 
  M <- sigY - U

  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL 
  
  ite <- 100
  epsilon <- 1e-5
  C1 <- diag(n)
  l2 <- rep(0, ite)
  Xn <- X
  Yn <- Y
  
	if (u == 0) {
	  for (i in 1:ite) {
	    m <- env(Xn, Yn, u, asy = F)
	    l2[i] <- m$loglik
	    resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
	    C1 <- diag(resid %*% solve(m$Sigma) %*% t(resid)) / r
	    Xn <- diag(1 / sqrt(C1)) %*% X
	    Yn <- diag(1 / sqrt(C1)) %*% Y
	    if (i > 1) {
	      if (abs(l2[i] - l2[i-1]) < epsilon * abs(l2[i-1])) break			
	    }
	  }	
	  Gammahat <- m$Gamma
	  Gamma0hat <- m$Gamma0
		etahat <- NULL
		Omegahat <- NULL
		Omega0hat <- m$Sigma
		muhat <- colMeans(Y)
		betahat <- matrix(0, r, p)
		Sigmahat <- m$Sigma
		loglik <- -n / 2 * sum(log(eigen(Sigmahat)$values)) - r * sum(log(C1)) / 2
		if (asy == T) ratio <- matrix(1, r, p)
		
	} else if (u == r) {
	
	  for (i in 1:ite) {
	    m <- env(Xn, Yn, u, asy = F)
	    l2[i] <- m$loglik
	    resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
	    C1 <- diag(resid %*% solve(m$Sigma) %*% t(resid)) / r
	    Xn <- diag(1 / sqrt(C1)) %*% X
	    Yn <- diag(1 / sqrt(C1)) %*% Y
	    if (i > 1) {
	      if (abs(l2[i] - l2[i-1]) < epsilon * abs(l2[i-1])) break			
	    }
	  }
	  Gammahat <- m$Gamma
	  Gamma0hat <- m$Gamma0
		etahat <- m$beta
		Omegahat <- diag(r)
		Omega0hat <- NULL
		muhat <- colMeans(Y) - m$beta %*% colMeans(X)
		betahat <- m$beta
		Sigmahat <- m$Sigma
		loglik <- -n / 2 * sum(log(eigen(Sigmahat)$values)) - r * sum(log(C1)) / 2
		if (asy == T) {
		  sigXtilde <- 0.25 * sigX
		  invsigXtilde <- chol2inv(chol(sigXtilde))
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigmahat)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- matrix(1, r, p)
		}
	} else {
		
		for (i in 1:ite) {
			m <- env(Xn, Yn, u, asy = F)
			l2[i] <- m$loglik
			resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
			C1 <- diag(resid %*% solve(m$Sigma) %*% t(resid)) / r
			Xn <- diag(1 / sqrt(C1)) %*% X
			Yn <- diag(1 / sqrt(C1)) %*% Y
			if (i > 1) {
				if (abs(l2[i] - l2[i-1]) < epsilon * abs(l2[i-1])) break			
			}
		}
		Gammahat <- m$Gamma
		Gamma0hat <- m$Gamma0
		muhat <- m$mu
		betahat <- m$beta
		etahat <- crossprod(Gammahat, betahat)
		Sigmahat <- m$Sigma
		Omegahat <- m$Omega
		Omega0hat <- m$Omega0
		loglik <- -n / 2 * sum(log(eigen(Sigmahat)$values)) - r * sum(log(C1)) / 2
		if (asy == T) {
		  sigXtilde <- 0.25 * sigX
		  invsigXtilde <- chol2inv(chol(sigXtilde))
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigmahat)
		  asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  invOmega0hat <- chol2inv(chol(Omega0hat))
		  invOmegahat <- chol2inv(chol(Omegahat))
	    M <- 0.25 * mean(C1)
		  temp <- kronecker(etahat %*% tcrossprod(sigXtilde, etahat), invOmega0hat) + M * (kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u)))
		  temp2 <- kronecker(t(etahat), Gamma0hat)
		  Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigma1) + 0.25 * temp2 %*% chol2inv(chol(temp)) %*% t(temp2)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- asyFm / asySE
		  }   
		}
 
   return(list(beta = betahat, Sigma = Sigmahat, Gamma = Gammahat, Gamma0 = Gamma0hat, eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, mu = muhat, loglik = loglik, covMatrix = covMatrix, asySE = asySE, ratio = ratio, n = n, C1 = C1))
}
	




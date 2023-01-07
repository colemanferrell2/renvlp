env.tcond <- function(X, Y, u, df, asy = TRUE) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)

  if (a[1] != nrow(X)) stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) stop("u must be an integer between 0 and r.")
  if(sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  
  Xc <- scale(X, center = T, scale = F)
  Yc <- scale(Y, center = T, scale = F)
  
  sigY <- stats::cov(Y) * (n - 1) / n
  sigYX <- stats::cov(Y, X) * (n - 1) / n
  sigX <- stats::cov(X) * (n - 1) / n
  invsigX <- chol2inv(chol(sigX))
  betaOLS <- sigYX %*% invsigX
  
  U <- tcrossprod(betaOLS, sigYX) 
  M <- sigY - U
  
  ite <- 100
  epsilon <- 1e-5
  C <- diag(n)
  l2 <- rep(0, ite)
  Xn <- X
  Yn <- Y
  
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL 
  
	if (u == 0) {
	
	  for (i in 1:ite) {
	    m <- env(Xn, Yn, u, asy = FALSE)
	    l2[i] <- m$loglik
	    resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
	    C1 <- diag(resid %*% solve(m$Sigma) %*% t(resid) + df) / (r + df)
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
	  alphahat <- colMeans(Y)
	  betahat <- matrix(0, r, p)
	  Sigmahat <- m$Sigma
		loglik <- -(r+df)/2 * sum(log(diag(Yc %*% solve(Sigmahat) %*% t(Yc)) / C1 / df + 1)) - 0.5 * n * sum(log(eigen(Sigmahat)$values)) - r * sum(log(C1)) / 2
		if (asy == T) ratio <- matrix(1, r, p)

		tmp <- Yc - Xc %*% t(m$beta)
		C1 <- (diag(tmp %*% solve(Sigmahat) %*% t(tmp)) + df) / (r + df) 
		loglik <- -(r+df)/2 * sum(log(diag(tmp %*% solve(m$Sigma) %*% t(tmp)) / C1 / df + 1)) - 0.5 * n * sum(log(eigen(m$Sigma)$values)) - r * sum(log(C1)) / 2
		
	} else {
		
		for (i in 1:ite) {
			m <- env(Xn, Yn, u, asy = FALSE)
			l2[i] <- m$loglik
			resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
			C1 <- (diag(resid %*% solve(m$Sigma) %*% t(resid)) + df) / (r + df) 
			Xn <- diag(1 / sqrt(C1)) %*% X
			Yn <- diag(1 / sqrt(C1)) %*% Y
			if (i > 1) {
				if (abs(l2[i] - l2[i-1]) < epsilon * abs(l2[i-1])) break			
			}
		}
		Gammahat <- m$Gamma
		Gamma0hat <- m$Gamma0
		alphahat <- m$mu
		betahat <- m$beta
		etahat <- crossprod(Gammahat, betahat)
		Sigmahat <- m$Sigma
		Omegahat <- m$Omega
		Omega0hat <- m$Omega0
		tmp <- Yc - Xc %*% t(m$beta)
		C1 <- (diag(tmp %*% solve(Sigmahat) %*% t(tmp)) + df) / (r + df) 
		loglik <- -(r+df)/2 * sum(log(diag(tmp %*% solve(m$Sigma) %*% t(tmp)) / C1 / df + 1)) - 0.5 * n * sum(log(eigen(m$Sigma)$values)) - r * sum(log(C1)) / 2
		if (asy == T) {
		  M <- (df + r)^2 / (4 * (r + df + 2) * (r + df))
		  Nx <- M
		  sigXtilde <- Nx * sigX
		  invsigXtilde <- chol2inv(chol(sigXtilde))
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigmahat)
		  if (u == r) {
		    asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		    ratio <- matrix(1, r, p)
		  } else {
		    asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
		    invOmega0hat <- chol2inv(chol(Omega0hat))
		    invOmegahat <- chol2inv(chol(Omegahat))
		    temp <- kronecker(etahat %*% tcrossprod(sigXtilde, etahat), invOmega0hat) + M * (kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u)))
		    temp2 <- kronecker(t(etahat), Gamma0hat)
		    Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
		    covMatrix <- 0.25 * kronecker(invsigXtilde, Sigma1) + 0.25 * temp2 %*% chol2inv(chol(temp)) %*% t(temp2)
		    asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		    ratio <- asyFm / asySE
		  }
		}
	}
 
   return(list(beta = betahat, Sigma = Sigmahat, Gamma = Gammahat, Gamma0 = Gamma0hat, eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, mu = alphahat, loglik = loglik, covMatrix = covMatrix, asySE = asySE, ratio = ratio, n = n))
}
	



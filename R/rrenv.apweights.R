rrenv.apweights <- function(X, Y, u, d, asy = TRUE) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  if (a[1] != nrow(X)) stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) stop("u must be an integer between 0 and r.")
  if(sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (d > u) stop("d must be no bigger than u.")
  
  
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
	    m <- env(Xn, Yn, u)
	    l2[i] <- m$loglik
	    resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
	    C1 <- diag(resid %*% solve(m$Sigma) %*% t(resid)) / r
	    Xn <- diag(1 / sqrt(C1)) %*% X
	    Yn <- diag(1 / sqrt(C1)) %*% Y
	    if (i > 1) {
	      if (abs(l2[i] - l2[i-1]) < epsilon * abs(l2[i-1])) break			
	    }
	  }
		Gammahat <- m$Gammahat
		Gamma0hat <- m$Gamma0hat
		etahat <- NULL
		Omegahat <- NULL
		Bhat <- NULL
		Omega0hat <- sigY
		alphahat <- colMeans(Y)
		betahat <- matrix(0, r, p)
		Sigmahat <- m$Sigma
		loglik <- -n / 2 * sum(log(eigen(Sigmahat)$values)) - r * sum(log(C1)) / 2
		if (asy == T) ratio <- matrix(1, r, p)
		
	} else if (u == r & d == r) {
	  for (i in 1:ite) {
	    m <- env(Xn, Yn, u)
	    l2[i] <- m$loglik
	    resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
	    C1 <- diag(resid %*% solve(m$Sigma) %*% t(resid)) / r
	    Xn <- diag(1 / sqrt(C1)) %*% X
	    Yn <- diag(1 / sqrt(C1)) %*% Y
	    if (i > 1) {
	      if (abs(l2[i] - l2[i-1]) < epsilon * abs(l2[i-1])) break			
	    }
	  }
		Gammahat <- m$Gammahat
		Gamma0hat <- m$Gamma0hat
		etahat <- m$beta
		Bhat <- diag(r)
		Omegahat <- diag(r)
		Omega0hat <- NULL
		alphahat <- colMeans(Y) - betaOLS %*% colMeans(X)
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
	} else if (u == r & d < r) {
	  for (i in 1:ite) {
	    m <- env(Xn, Yn, u)
        sigX <- stats::cov(Xn) * (n - 1) / n
        sigY <- stats::cov(Yn) * (n - 1) / n
        sigYX <- stats::cov(Yn, Xn) * (n - 1) / n
	    eigSx <- eigen(sigX)
	    Sxnhalf <- sweep(eigSx$vectors, MARGIN = 2, 1 / sqrt(eigSx$values), '*') %*% t(eigSx$vectors) 
	    eigsigY <- eigen(sigY)
	    sigYnhalf <- sweep(eigsigY$vectors, MARGIN = 2, 1 / sqrt(eigsigY$values), '*') %*% t(eigsigY$vectors)
	    sigYhalf <- sweep(eigsigY$vectors, MARGIN = 2, sqrt(eigsigY$values), '*') %*% t(eigsigY$vectors)
	    CYX <- sigYnhalf %*% sigYX %*% Sxnhalf
	    svdC <- svd(CYX)
	    if (d == 1) {
	      svdCd <- svdC$d[1:d] * as.matrix(svdC$u[, 1:d]) %*% t(as.matrix(svdC$v[, 1:d]))
	    } else {
	      svdCd <- svdC$u[, 1:d] %*% diag(svdC$d[1:d]) %*% t(svdC$v[, 1:d])
	    }
	    A <- sigYhalf %*% sweep(svdC$u[, 1:d], MARGIN = 2, svdC$d[1:d], '*')
	    B <- t(svdC$v[, 1:d]) %*% Sxnhalf
	    Bhat <- B
	    betahat <- A %*% B
	    etahat <- A
#	    betahat <- sigYhalf %*% svdCd %*% Sxnhalf
	    alphahat <- colMeans(Y) - betahat %*% colMeans(X)
	    Sigmahat <- sigYhalf %*% (diag(r) - svdCd %*% t(svdCd)) %*% sigYhalf
	    Omegahat <- Sigmahat
	    Omega0hat <- NULL
	    resid <- Yn - as.matrix(rep(1, n)) %*% t(alphahat) - Xn %*% t(betahat)
	    C1 <- diag(resid %*% solve(Sigmahat) %*% t(resid)) / r
	    l2[i] <- - n * r / 2 * (log(2 * pi) + 1) - n / 2 * sum(log(eigen(Sigmahat)$values)) 
	    Xn <- diag(1 / sqrt(C1)) %*% X
	    Yn <- diag(1 / sqrt(C1)) %*% Y
	    if (i > 1) {
	      if (abs(l2[i] - l2[i-1]) < epsilon * abs(l2[i-1])) break			
	    }
	  }
	
		Gammahat <- m$Gammahat
		Gamma0hat <- m$Gamma0hat
		loglik <- -n / 2 * sum(log(eigen(Sigmahat)$values)) - r * sum(log(C1)) / 2
		if (asy == T) {
		  sigXtilde <- 0.25 * sigX
		  invsigXtilde <- chol2inv(chol(sigXtilde))
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigmahat)
		  asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  invSigmahat <- chol2inv(chol(Sigmahat))
		  MA <- A %*% solve(crossprod(A, invSigmahat) %*% A) %*% t(A)
		  MB <- t(B) %*% solve(B %*% sigXtilde %*% t(B)) %*% B
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigmahat) - 0.25 * kronecker(invsigXtilde - MB, Sigmahat - MA)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- asyFm / asySE
		}  
		
	} else if (d == u | (d == p & p < r)) {
	
		tmp <- env.apweights(X, Y, u)
		Gammahat <- tmp$Gamma
		Gamma0hat <- tmp$Gamma0	
		betahat <- tmp$beta
		etahat <- tmp$eta
		Bhat <- diag(p)
		alphahat <- tmp$alpha
		Omegahat <- tmp$Omega
		Omega0hat <- tmp$Omega0
		Sigmahat <- tmp$Sigma
		loglik <- tmp$loglik
		C1 <- tmp$C1
		covMatrix <- tmp$covMatrix
		asySE <- tmp$asySE
		ratio <- tmp$ratio
	
	} else {

		for (i in 1:ite) {
			m <- rrenv(Xn, Yn, u, d)
			l2[i] <- m$loglik
			resid <- Yn - as.matrix(rep(1, n)) %*% t(m$mu) - Xn %*% t(m$beta)
			C1 <- diag(resid %*% chol2inv(chol(m$Sigma)) %*% t(resid)) / r
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
		etahat <- m$eta
		Bhat <- m$B
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
		  sigXtilde <- 0.25 * sigX
		  invsigXtilde <- chol2inv(chol(sigXtilde))
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigmahat)
		  asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  
		  invOmega0hat <- chol2inv(chol(Omega0hat))
		  invOmegahat <- chol2inv(chol(Omegahat))
		  tmp1 <- etahat %*% Bhat
		  tmp2 <- kronecker(t(tmp1), Gamma0hat)
		  tmp3 <- kronecker(tmp1 %*% sigXtilde %*% t(tmp1), invOmega0hat) + M * (kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u)))
		  tmp4 <- Gammahat %*% etahat %*% solve(t(etahat) %*% invOmegahat %*% etahat) %*% t(etahat) %*% t(Gammahat)
		  MB <- t(Bhat) %*% solve(Bhat %*% sigXtilde %*% t(Bhat)) %*% Bhat
		  covMatrix <- 0.25 * kronecker(invsigXtilde, Sigmahat) - 0.25 * kronecker(invsigXtilde - MB, Sigmahat - tmp4) - 0.25 * kronecker(MB, Gamma0hat %*% Omega0hat %*% t(Gamma0hat)) + 0.25 * tmp2 %*% solve(tmp3) %*% t(tmp2)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- asyFm / asySE
		}   
	}
 
   return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = alphahat, beta = betahat, Sigma = Sigmahat, eta = etahat, B = Bhat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik, covMatrix = covMatrix, asySE = asySE, ratio = ratio, n = n, C1 = C1))
}
	



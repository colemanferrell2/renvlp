rrenv <- function(X, Y, u, d, asy = TRUE) {

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
  M <- as.matrix(Matrix::nearPD(M)$mat)

  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL 
  
	if (u == 0) {
	
	  tmp <- envMU(M, U, u)
	  Gammahat <- tmp$Gammahat
	  Gamma0hat <- tmp$Gamma0hat
		etahat <- NULL
		Bhat <- NULL
		Omegahat <- NULL
		Omega0hat <- sigY
		alphahat <- colMeans(Y)
		betahat <- matrix(0, r, p)
		Sigmahat <- sigY
		if (asy == T) ratio <- matrix(1, r, p)
		
	} else if (u == r & d == r) {
	
	  tmp <- envMU(M, U, u)
	  Gammahat <- tmp$Gammahat
	  Gamma0hat <- tmp$Gamma0hat
		etahat <- betaOLS
		Bhat <- diag(r)
		Omegahat <- diag(r)
		Omega0hat <- NULL
		alphahat <- colMeans(Y) - betaOLS %*% colMeans(X)
		betahat <- betaOLS
		Sigmahat <- M
		if (asy == T) {
		  covMatrix <- kronecker(invsigX, M)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- matrix(1, r, p)
		}		
		
	} else if (u == r & d < r) {
	
	  tmp <- envMU(M, U, u)
	  Gammahat <- tmp$Gammahat
	  Gamma0hat <- tmp$Gamma0hat
		eigSx <- eigen(sigX)
		Sxnhalf <- sweep(eigSx$vectors, MARGIN = 2, 1 / sqrt(eigSx$values), '*') %*% t(eigSx$vectors) 
		eigsigY <- eigen(sigY)
		sigYnhalf <- sweep(eigsigY$vectors, MARGIN = 2, 1 / sqrt(eigsigY$values), '*') %*% t(eigsigY$vectors)
		sigYhalf <- sweep(eigsigY$vectors, MARGIN = 2, sqrt(eigsigY$values), '*') %*% t(eigsigY$vectors)
		CYX <- sigYnhalf %*% sigYX %*% Sxnhalf
		svdC <- svd(CYX)
		svdCd <- sweep(svdC$u[, 1:d], MARGIN = 2, svdC$d[1:d], '*') %*% t(svdC$v[, 1:d])
		A <- sigYhalf %*% sweep(svdC$u[, 1:d], MARGIN = 2, svdC$d[1:d], '*')
		B <- t(svdC$v[, 1:d]) %*% Sxnhalf
		Bhat <- B
		betahat <- A %*% B
#		betahat <- sigYhalf %*% svdCd %*% Sxnhalf
		etahat <- A
		alphahat <- colMeans(Y) - betahat %*% colMeans(X)
		Sigmahat <- sigYhalf %*% (diag(r) - svdCd %*% t(svdCd)) %*% sigYhalf
		Omegahat <- Sigmahat
		Omega0hat <- NULL
		if (asy == T) {
		  covMatrix <- kronecker(invsigX, Sigmahat)
		  asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  invSigmahat <- chol2inv(chol(Sigmahat))
		  MA <- A %*% solve(crossprod(A, invSigmahat) %*% A) %*% t(A)
		  MB <- t(B) %*% solve(B %*% sigX %*% t(B)) %*% B
		  covMatrix <- kronecker(invsigX, Sigmahat) - kronecker(invsigX - MB, Sigmahat - MA)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- asyFm / asySE
		}  

	} else if (d == u | (d == p & p < r)) {
	
	  tmp <- envMU(M, U, u)
	  Gammahat <- tmp$Gammahat
	  Gamma0hat <- tmp$Gamma0hat
		etahat <- crossprod(Gammahat, betaOLS)
		Bhat <- diag(p)
		betahat <- Gammahat %*% etahat
		alphahat <- colMeans(Y) - betahat %*% colMeans(X)
		Omegahat <- crossprod(Gammahat, M) %*% Gammahat
		Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
		Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
		Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)	
		if (asy == T) {
		  covMatrix <- kronecker(invsigX, Sigmahat)
		  asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  invOmega0hat <- chol2inv(chol(Omega0hat))
		  invOmegahat <- chol2inv(chol(Omegahat))
		  tmp2 <- kronecker(t(etahat), Gamma0hat)
		  tmp3 <- kronecker(etahat %*% sigX %*% t(etahat), invOmega0hat) + kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u))
		  covMatrix <- kronecker(invsigX, Gammahat %*% Omegahat %*% t(Gammahat)) + tmp2 %*% solve(tmp3) %*% t(tmp2)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- asyFm / asySE
		}  	
	} else {
		tmp <- rrenvMU(M, U, u, d)
		Gammahat <- tmp$Gammahat
		Gamma0hat <- tmp$Gamma0hat
		eigSx <- eigen(sigX)
		Sxnhalf <- sweep(eigSx$vectors, MARGIN = 2, 1 / sqrt(eigSx$values), '*') %*% t(eigSx$vectors) 
		SgamY <- crossprod(Gammahat, sigY) %*% Gammahat
		eigSgamY <- eigen(SgamY)
		SgamYnhalf <- sweep(eigSgamY$vectors, MARGIN = 2, 1 / sqrt(eigSgamY$values), '*') %*% t(eigSgamY$vectors)
		SgamYhalf <- sweep(eigSgamY$vectors, MARGIN = 2, sqrt(eigSgamY$values), '*') %*% t(eigSgamY$vectors)
		CgamYX <- SgamYnhalf %*% t(Gammahat) %*% sigYX %*% Sxnhalf
		svdC <- svd(CgamYX)
		svdCd <- sweep(as.matrix(svdC$u[, 1:d]), MARGIN = 2, svdC$d[1:d], '*') %*% t(svdC$v[, 1:d])
		etahat <- SgamYhalf %*% sweep(as.matrix(svdC$u[, 1:d]), MARGIN = 2, svdC$d[1:d], '*')
		Bhat <- t(svdC$v[, 1:d]) %*% Sxnhalf
		betahat <- Gammahat %*% SgamYhalf %*% svdCd %*% Sxnhalf
		alphahat <- colMeans(Y) - betahat %*% colMeans(X)
		Omegahat <- SgamYhalf %*% (diag(u) - svdCd %*% t(svdCd)) %*% SgamYhalf
		Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
		Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
		Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
		if (asy == T) {
		  covMatrix <- kronecker(invsigX, Sigmahat)
		  asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  invOmega0hat <- chol2inv(chol(Omega0hat))
		  invOmegahat <- chol2inv(chol(Omegahat))
		  tmp1 <- etahat %*% Bhat
		  tmp2 <- kronecker(t(tmp1), Gamma0hat)
		  tmp3 <- kronecker(tmp1 %*% sigX %*% t(tmp1), invOmega0hat) + kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u))
		  tmp4 <- Gammahat %*% etahat %*% solve(t(etahat) %*% invOmegahat %*% etahat) %*% t(etahat) %*% t(Gammahat)
		  MB <- t(Bhat) %*% solve(Bhat %*% sigX %*% t(Bhat)) %*% Bhat
		  covMatrix <- kronecker(invsigX, Sigmahat) - kronecker(invsigX - MB, Sigmahat - tmp4) - kronecker(MB, Gamma0hat %*% Omega0hat %*% t(Gamma0hat)) + tmp2 %*% solve(tmp3) %*% t(tmp2)
		  asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
		  ratio <- asyFm / asySE
		}  
	}
	
	t1 <- eigen(Sigmahat)
	invSigmahat <- sweep(t1$vectors, MARGIN = 2, 1 / t1$values, '*') %*% t(t1$vectors)
	loglik <- -n / 2 * sum(log(t1$values)) - n / 2 * sum(diag(invSigmahat %*% (M + (betaOLS - betahat) %*% sigX%*% t(betaOLS - betahat))))
 
	return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = alphahat, beta = betahat, Sigma = Sigmahat, eta = etahat, B = Bhat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik, covMatrix = covMatrix, asySE = asySE, ratio = ratio, n = n))
}
	



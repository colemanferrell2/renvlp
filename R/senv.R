senv <- function(X, Y, u, asy = TRUE, init = NULL){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  if (a[1] != nrow(X)) 
    stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) 
    stop("u must be an integer between 0 and r.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
    if (nrow(init) != r || ncol(init) != u) stop("The dimension of init is wrong.")
  }
  
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigX <- chol2inv(chol(sigX))
  betaOLS <- sigYX %*% invsigX
  U <- tcrossprod(betaOLS, sigYX)
  M <- sigY - U
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL
  if (u == 0) {
    Gammahat <- NULL
    Gamma0hat <- diag(r)
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigY
    Lambdahat <- diag(r)
    muhat <- colMeans(Y)
    betahat <- matrix(0, r, p)
    Sigmahat <- sigY
    tmp.MU <- eigen(sigY)
    objfun <- sum(log(tmp.MU$values))
    loglik <- -n * r/2 * (log(2 * pi) + 1) - n/2 * objfun
    if (asy == T) 
      ratio <- matrix(1, r, p)
  }
  else if (u >= (r - (r-1)/p)) {
    Gammahat <- diag(r)
    Gamma0hat <- NULL
    etahat <- betaOLS
    Omegahat <- M
    Omega0hat <- NULL
    Lambdahat <- diag(r)
    muhat <- colMeans(Y) - betaOLS %*% colMeans(X)
    betahat <- betaOLS
    Sigmahat <- M
    tmp.M <- eigen(M)
    objfun <- sum(log(tmp.M$values))
    loglik <- -n * r/2 * (log(2 * pi) + 1) - n/2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(invsigX, M)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- matrix(1, r, p)
    }
  }
  else {
    tmp <- senvMU(X, Y, u, initial = init)
    
    
    Gammahat <- tmp$Gammahat
    Gamma0hat <- tmp$Gamma0hat
    Lambdahat <- tmp$Lambda
    d1 <- diag(Lambdahat)
    invLambda <- diag(1 / d1)
    etahat <- crossprod(Gammahat, invLambda) %*% betaOLS
    betahat <- Lambdahat %*% Gammahat %*% etahat
    muhat <- colMeans(Y) - betahat %*% colMeans(X)
    M1 <- invLambda %*% M %*% invLambda
    Omegahat <- crossprod(Gammahat, M1) %*% Gammahat
    M2 <- invLambda %*% sigY %*% invLambda
    Omega0hat <- crossprod(Gamma0hat, M2) %*% Gamma0hat
    E1 <- Lambdahat %*% Gammahat
    Sigma1 <- E1 %*% tcrossprod(Omegahat, E1)
    E2 <- Lambdahat %*% Gamma0hat
    Sigmahat <- Sigma1 + E2 %*% tcrossprod(Omega0hat, E2)
    tmp.MU <- eigen(sigY)
    objfun <- tmp$objfun + sum(log(tmp.MU$values))
    loglik <- -n * r/2 * (log(2 * pi) + 1) - n/2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(invsigX, M)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
      SigmaL <- Gammahat %*% tcrossprod(Omegahat, 
                Gammahat) + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
      invsig <- chol2inv(chol(Sigmahat))
      Cr <- contr(r)
      Er <- expan(r)
      Eu <- expan(u)
      Eru <- expan(r - u) 
      g11 <- kronecker(diag(p), Gammahat)
      g12 <- kronecker(t(etahat), diag(r))
      g221 <- kronecker(Gammahat %*% Omegahat, diag(r))
      g222 <- kronecker(Gammahat, Gamma0hat %*% 
                          tcrossprod(Omega0hat, Gamma0hat))
      g22 <- 2 * Cr %*% (g221 - g222)
      g23 <- Cr %*% kronecker(Gammahat, Gammahat) %*% Eu
      g24 <- Cr %*% kronecker(Gamma0hat, Gamma0hat) %*% Eru
      G_o <- cbind(rbind(g11, matrix(0, nrow = nrow(g22), ncol = ncol(g11))),
                   rbind(g12, g22),
                   rbind(matrix(0, nrow = nrow(g11), ncol = ncol(g23)), g23),
                   rbind(matrix(0, nrow = nrow(g11), ncol = ncol(g24)), g24))
      h_o1 <- kronecker(Gammahat %*% etahat, diag(r))
      h_o2 <- 2 * tcrossprod(kronecker(SigmaL, diag(r)), Cr) 
      th_o <- cbind(h_o1, h_o2)
      h_o <- t(th_o)
      d11 <- kronecker(diag(p), Lambdahat)
      d22 <- Cr %*% kronecker(Lambdahat, Lambdahat) %*% Er
      DL <- cbind(rbind(d11, matrix(0, nrow = nrow(d22), ncol = ncol(d11))), 
                  rbind(matrix(0, nrow = nrow(d11), ncol = ncol(d22)), d22))
      diagvec <- function(d){
        E <- NULL
        for (i in 1 : (d - 1)){
          e <- rep(0, d)
          e[i + 1] <- 1
          t <- kronecker(e, e)
          E <- cbind(E, t)
        }
        E <- matrix(E, nrow = d^2)
      }
      L <- diagvec(r)
      h1 <- DL %*% h_o %*% kronecker(diag(r), invLambda) %*% L
      h2 <- DL %*% G_o
      H <- cbind(h1, h2)
      j11 <- kronecker(sigX, invsig)
      j221 <- kronecker(invsig, invsig)
      j22 <- 1/2 * crossprod(Er, j221) %*% Er
      J <- cbind(rbind(j11, matrix(0, nrow = nrow(j22), ncol = ncol(j11))),
                 rbind(matrix(0, nrow = nrow(j11), ncol = ncol(j22)), j22))
      M1 <- crossprod(H, J) %*% H
      invM1 <- ginv(M1)
      V1 <- tcrossprod(invM1, H)
      V <- H %*% V1
      covMatrix <- V[1:(p*r), 1:(p*r)]
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- asyFm/asySE
    }
  }  
  return(list(beta = betahat, Sigma = Sigmahat, Lambda = Lambdahat,
              Gamma = Gammahat, Gamma0 = Gamma0hat, eta = etahat, 
              Omega = Omegahat, Omega0 = Omega0hat, mu = muhat,
              loglik = loglik, covMatrix = covMatrix, 
              asySE = asySE, ratio = ratio, n = n))
}

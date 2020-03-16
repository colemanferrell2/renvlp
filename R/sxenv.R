sxenv <- function(X, Y, u, R, asy = TRUE, init = NULL){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  if (a[1] != nrow(X)) 
    stop("X and Y should have the same number of observations.")
  if (u > p || u < 0) 
    stop("u must be an integer between 0 and p.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  if (!is.null(init)) {
    if (nrow(init) != p || ncol(init) != u) stop("The dimension of init is wrong.")
  }
  
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigY <- chol2inv(chol(sigY))
  invsigX <- chol2inv(chol(sigX))
  eig.sigX <- eigen(sigX)
  eig.sigY <- eigen(sigY)
  U <- crossprod(sigYX, invsigY) %*% sigYX
  M <- sigX - U
  q <- length(R)
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL
  if (u == 0) {
    Gammahat <- NULL
    Gamma0hat <- diag(p)
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigX
    Lambdahat <- diag(p)
    muYhat <- colMeans(Y)
    muXhat <- colMeans(X)
    betahat <- matrix(0, p, r)
    SigmaXhat <- sigX
    SigmaYcXhat <- sigY
    objfun <- sum(log(eig.sigX$values)) + sum(log(eig.sigY$values)) + p
    loglik <- -n * (p + r)/2 * log(2 * pi) - n * r/2 - n/2 * objfun
    if (asy == T) 
      ratio <- matrix(1, p, r)
  } else if (u >= (p - (q - 1) / r)) {
    Gammahat <- diag(p)
    Gamma0hat <- NULL
    invsigX <- chol2inv(chol(sigX))
    betaOLS <- tcrossprod(invsigX, sigYX)
    etahat <- betaOLS
    Lambdahat <- diag(p)
    Omegahat <- sigX
    Omega0hat <- NULL
    muYhat <- colMeans(Y)
    muXhat <- colMeans(X)
    betahat <- betaOLS
    SigmaXhat <- sigX
    SigmaYcXhat <- sigY - sigYX %*% betaOLS
    eig.sigXhat <- eigen(SigmaXhat)
    eig.sigYcXhat <- eigen(SigmaYcXhat)
    objfun <- (sum(log(eig.sigXhat$values)) + sum(log(eig.sigYcXhat$values)) 
               + p)
    loglik <- -n * (r + p)/2 * log(2 * pi) - n * r/2 - n/2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(SigmaYcXhat, invsigX)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- matrix(1, p, r)
    }
  } else {
    tmp1 <- sxenvMU(X, Y, u, R, initial = init)
    
    Gammahat <- tmp1$Gammahat
    Gamma0hat <- tmp1$Gamma0hat
    Lambdahat <- tmp1$Lambdahat
    d1 <- diag(Lambdahat)
    invLambdahat <- diag(1 / d1)
    E1 <- crossprod(Gammahat, invLambdahat)
    Omegahat <- E1 %*% tcrossprod(sigX, E1)
    E2 <- crossprod(Gamma0hat, invLambdahat)
    Omega0hat <- E2 %*% tcrossprod(sigX, E2)
    invOmegahat <- chol2inv(chol(Omegahat))
    etahat <- invOmegahat %*% tcrossprod(E1, sigYX) 
    betahat <- invLambdahat %*% Gammahat %*% etahat
    muYhat <- colMeans(Y)
    muXhat <- colMeans(X)
    E3 <- Lambdahat %*% Gammahat
    sig1 <- E3 %*% tcrossprod(Omegahat, E3)
    E4 <- Lambdahat %*% Gamma0hat
    sig2 <- E4 %*% tcrossprod(Omega0hat, E4)
    SigmaXhat <- sig1 + sig2
    Yc <- Y - tcrossprod(rep(1, nrow(Y)), colMeans(Y))
    Xc <- X - tcrossprod(rep(1, nrow(X)), colMeans(X))
    E5 <- Yc - Xc %*% betahat
    SigmaYcXhat <- crossprod(E5) / n
    invsigYcX <- chol2inv(chol(SigmaYcXhat))
    eig.sigXhat <- eigen(SigmaXhat)
    eig.sigYcXhat <- eigen(SigmaYcXhat)
    invsigXhat <- chol2inv(chol(SigmaXhat))
    s1 <- invsigXhat %*% sigX
    eig.invX <- eigen(s1)
    objfun <- (sum(log(eig.sigXhat$values)) + sum(log(eig.sigYcXhat$values)) 
               + sum(eig.invX$values))
    loglik <- -n * (r + p)/2 * log(2 * pi) - n * r/2 - n/2 * objfun
    if (asy == T) {
      U1 <- sigYX %*% invsigX
      U2 <- tcrossprod(U1, sigYX)
      M2 <- sigY - U2
      covMatrix <- kronecker(M2, invsigX)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = p)
      Ep <- expan(p)
      Eu <- expan(u)
      Epu <- expan((p - u))
      Cp <- contr(p)
      j11 <- kronecker(invsigYcX, SigmaXhat)
      j222 <- kronecker(invsigXhat, invsigXhat)
      j22 <- crossprod(Ep, j222) %*% Ep / 2
      J <- cbind(rbind(j11, matrix(0, nrow = nrow(j22), ncol = ncol(j11))), 
                  rbind(matrix(0, nrow = nrow(j11), ncol = ncol(j22)), j22))
      diagvec <- function(R){
        E <- NULL
        p <- sum(R)
        for(i in 1 : (length(R) - 1)){
          e <- rep(0, p)
          s <- sum(R[1 : i])
          e[s + 1] <- 1
          t <- kronecker(e, e)
          E <- cbind(E, t)
        }
        E <- matrix(E, nrow = p^2)
      }
      L <- diagvec(R)
      SigmaL <- Gammahat %*% tcrossprod(Omegahat, 
                Gammahat) + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
      h1 <- crossprod(Gammahat, invLambdahat)
      h2 <- crossprod(etahat, h1)
      h11 <- - kronecker(h2, invLambdahat) %*% L
      h3 <- kronecker(Lambdahat %*% SigmaL, diag(p))
      h21 <- 2 * Cp %*% h3 %*% L
      h12 <- kronecker(diag(r), invLambdahat %*% Gammahat)
      h13 <- kronecker(t(etahat), invLambdahat)
      h4 <- Lambdahat %*% Gammahat %*% Omegahat
      h5 <- kronecker(h4, Lambdahat)
      h6 <- Lambdahat %*% Gammahat
      h7 <- Lambdahat %*% Gamma0hat 
      h78 <- h7 %*% tcrossprod(Omega0hat, Gamma0hat)
      h8  <- kronecker(h6, h78)
      h23 <- 2 * Cp %*% (h5 - h8)
      h9 <- Lambdahat %*% Gammahat
      h24 <- Cp %*% kronecker(h9, h9) %*% Eu
      h10 <- Lambdahat %*% Gamma0hat
      h25 <- Cp %*% kronecker(h10, h10) %*% Epu
      H <- rbind(cbind(h11, h12, h13, matrix(0, nrow = nrow(h11), 
                ncol = ncol(h24)), matrix(0, nrow = nrow(h11), ncol = ncol(h25))),
                cbind(h21, matrix(0, nrow = nrow(h21), ncol = ncol(h12)),
                      h23, h24, h25))
      M1 <- crossprod(H, J) %*% H
      invM1 <- ginv(M1)
      V <- H %*% tcrossprod(invM1, H)
      covMatrix <- V[1:(p * r), 1:(p * r)]
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- asyFm/asySE
    }
  }
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat,
              muY = muYhat, muX = muXhat, Lambda = Lambdahat, 
              beta = as.matrix(betahat), SigmaX = SigmaXhat, eta = etahat, 
              Omega = Omegahat, Omega0 = Omega0hat, SigmaYcX = SigmaYcXhat, 
              loglik = loglik, n = n, covMatrix = covMatrix, asySE = asySE, 
              ratio = ratio))
}

stenv <- function (X, Y, q, u, asy = TRUE, Pinit = NULL, Ginit = NULL) {
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
  if (q > p || q < 0) 
    stop("q must be an integer between 0 and p.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigX <- chol2inv(chol(sigX))
  betaOLS <- sigYX %*% invsigX
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL
  if (u == 0 || q ==0) {
    Gammahat <- NULL
    Gamma0hat <- diag(r)
    Phihat <- NULL
    Phi0hat <- diag(p)
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigY
    Deltahat <- NULL
    Delta0hat <- sigX
    muhat <- colMeans(Y)
    betahat <- matrix(0, p, r)
    SigmaYcXhat <- sigY
    SigmaXhat <- sigX
    tmp1 <- eigen(sigX)
    tmp2 <- eigen(sigY)
    objfun <- sum(log(tmp1$values)) + sum(log(tmp2$values))
    loglik <- - n/2 * objfun
    if (asy == T) 
      ratio <- matrix(1, p, r)
  }
  else if (u == r && q == p) {
    M <- sigY - tcrossprod(betaOLS, sigYX)
    Gammahat <- diag(r)
    Gamma0hat <- NULL
    Phihat <- diag(p)
    Phi0hat <- NULL
    etahat <- betaOLS
    Omegahat <- M
    Omega0hat <- NULL
    Deltahat <- sigX
    Delta0hat <- NULL
    muhat <- colMeans(Y) - betaOLS %*% colMeans(X)
    betahat <- t(betaOLS)
    SigmaYcXhat <- M
    SigmaXhat <- sigX
    tmp.e1 <- eigen(Omegahat)
    tmp.e3 <- eigen(Deltahat)
    objfun <- sum(log(tmp.e1$values)) + sum(log(tmp.e3$values)) 
    loglik <- -n * (p + r) / 2 - n/2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(M, invsigX)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- matrix(1, p, r)
    }
  }
  else if (u == r) {
    M1 <- sigY - tcrossprod(betaOLS, sigYX)
    tmp.xenv <- xenv(X, Y, q, asy = TRUE)
    Gammahat <- diag(r)
    Gamma0hat <- NULL
    Phihat <- tmp.xenv$Gamma
    Phi0hat <- tmp.xenv$Gamma0
    E1 <-  crossprod(Phihat, sigX) %*% Phihat
    invE1 <- chol2inv(chol(E1))
    etahat <- invE1 %*% t(Phihat) %*% t(sigYX)
    Omegahat <- M1
    Omega0hat <- NULL
    Deltahat <- tmp.xenv$Omega
    Delta0hat <- tmp.xenv$Omega0
    muhat <- tmp.xenv$mu
    betahat <- Phihat %*% etahat
    SigmaYcXhat <- tmp.xenv$SigmaYcX
    SigmaXhat <- tmp.xenv$SigmaX
    invome <- chol2inv(chol(Omegahat))
    invdel <- chol2inv(chol(Deltahat))
    invdel0 <- chol2inv(chol(Delta0hat))
    tmp.e1 <- eigen(Omegahat)
    tmp.e2 <- eigen(Deltahat)
    tmp.e3 <- eigen(Delta0hat)
    t3 <- sigY %*% invome
    t5 <- crossprod(betahat, sigX) %*% betahat %*% invome
    t6 <- t(betahat) %*% t(sigYX) %*% invome
    obj1 <- sum(log(tmp.e1$values)) + sum(log(tmp.e2$values)) 
    obj2 <- sum(log(tmp.e3$values)) + p
    obj3 <- sum(diag(t3))
    obj4 <- sum(diag(t5)) - 2 * sum(diag(t6))
    objfun <- obj1 + obj2 + obj3 + obj4
    loglik <- - (n/2) * objfun
    if (asy == T) {
      covMatrix <- tmp.xenv$covMatrix
      asySE <- tmp.xenv$asySE
      ratio <- tmp.xenv$ratio
    }
  }
  else if (q == p) {
    tmp.env <- env(X, Y, u, asy = TRUE)
    Gammahat <- tmp.env$Gamma
    Gamma0hat <- tmp.env$Gamma0
    Phihat <- diag(p)
    Phi0hat <- NULL
    etahat <- invsigX %*% crossprod(sigYX, Gammahat)
    Omegahat <- tmp.env$Omega
    Omega0hat <- tmp.env$Omega0
    Deltahat <- sigX
    Delta0hat <- NULL
    betahat <- t(tmp.env$beta)
    muhat <- colMeans(Y) - crossprod(betahat, colMeans(X))
    SigmaYcXhat <- tmp.env$Sigma
    SigmaXhat <- sigX
    tmp.e1 <- eigen(Omegahat)
    tmp.e2 <- eigen(Omega0hat)
    tmp.e3 <- eigen(Deltahat)
    invome <- chol2inv(chol(Omegahat))
    invome0 <- chol2inv(chol(Omega0hat))
    invdel <- chol2inv(chol(Deltahat))
    t3 <- crossprod(Gammahat, sigY) %*% Gammahat %*% invome
    t5 <- crossprod(etahat, sigX) %*% etahat %*% invome
    t6 <- t(etahat) %*% crossprod(sigYX, Gammahat) %*% invome
    obj1 <- sum(log(tmp.e1$values)) + sum(log(tmp.e2$values)) 
    obj2 <- sum(log(tmp.e3$values)) + p
    obj3 <- sum(diag(t3)) + (r - u)
    obj4 <- sum(diag(t5)) - 2 * sum(diag(t6))
    objfun <- obj1 + obj2 + obj3 + obj4
    loglik <- - (n/2) * objfun
    if (asy == T) {
      M <- sigY - tcrossprod(betaOLS, sigYX)
      covMatrix <- kronecker(M, invsigX)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = p)
      m1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
      temp1 <- kronecker(m1, invsigX)
      m2 <- crossprod(etahat, sigX) %*% etahat
      m3 <- kronecker(invome0, m2) + kronecker(invome0, Omegahat)
      m4 <- kronecker(Omega0hat, invome) - 2 * kronecker(diag(r - u), diag(u))
      m5 <- m3 + m4
      invm5 <- chol2inv(chol(m5))
      m6 <- kronecker(Gamma0hat, etahat)
      temp2 <- m6 %*% tcrossprod(invm5, m6)
      covMatrix <- temp1 + temp2
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- asyFm/asySE
    }
  }
  else{

    tmp.stenv <- stenvMU(X, Y, q, u)

    if (!is.null(Ginit)) {
        if (nrow(Ginit) != r || ncol(Ginit) != u) stop("The initial value should have r rows and u columns.")
        tmp0 <- qr.Q(qr(Ginit), complete = TRUE)
        tmp.stenv$Gammahat <- as.matrix(tmp0[, 1:u])
        tmp.stenv$Gamma0hat <- as.matrix(tmp0[, (u+1):r])
    }
    
    if (!is.null(Pinit)) {
        if (nrow(Pinit) != p || ncol(Pinit) != q) stop("The initial value should have p rows and q columns.")
        tmp0 <- qr.Q(qr(Pinit), complete = TRUE)
        tmp.stenv$Phihat <- as.matrix(tmp0[, 1:q])
        tmp.stenv$Phi0hat <- as.matrix(tmp0[, (q+1):p])
    }
    
    Gammahat <- tmp.stenv$Gammahat
    Gamma0hat <- tmp.stenv$Gamma0hat
    Phihat <- tmp.stenv$Phihat
    Phi0hat <- tmp.stenv$Phi0hat
    Deltahat <- crossprod(Phihat, sigX) %*% Phihat
    Delta0hat <- crossprod(Phi0hat, sigX) %*% Phi0hat
    invdel <- chol2inv(chol(Deltahat))
    et1 <- t(Phihat) %*% crossprod(sigYX, Gammahat)
    etahat <- invdel %*% et1
    om1 <- sigYX %*% Phihat
    om2 <- om1 %*% tcrossprod(invdel, om1)
    om3 <- sigY - om2
    Omegahat <- crossprod(Gammahat, om3) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
    betahat <- Phihat %*% tcrossprod(etahat, Gammahat)
    muhat <- colMeans(Y) - crossprod(betahat, colMeans(X))
    s1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
    SigmaYcXhat <- s1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
    s2 <- Phihat %*% tcrossprod(Deltahat, Phihat)
    SigmaXhat <- s2 + Phi0hat %*% tcrossprod(Delta0hat, Phi0hat)
    tmp.e1 <- eigen(Omegahat)
    tmp.e2 <- eigen(Omega0hat)
    tmp.e3 <- eigen(Deltahat)
    tmp.e4 <- eigen(Delta0hat)
    invome <- chol2inv(chol(Omegahat))
    invome0 <- chol2inv(chol(Omega0hat))
    invdel0 <- chol2inv(chol(Delta0hat))
    t3 <- crossprod(Gammahat, sigY) %*% Gammahat %*% invome
    peta <- Phihat %*% etahat
    t5 <- crossprod(peta, sigX) %*% peta %*% invome
    t6 <- t(peta) %*% crossprod(sigYX, Gammahat) %*% invome
    obj1 <- sum(log(tmp.e1$values)) + sum(log(tmp.e2$values)) 
    obj2 <- sum(log(tmp.e3$values)) + sum(log(tmp.e4$values)) + p
    obj3 <- sum(diag(t3)) + (r - u)
    obj4 <- sum(diag(t5)) - 2 * sum(diag(t6))
    objfun <- obj1 + obj2 + obj3 + obj4
    loglik <- - (n/2) * objfun
    if (asy == T) {
      M <- sigY - tcrossprod(betaOLS, sigYX)
      covMatrix <- kronecker(M, invsigX)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = p)
      m1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
      m2 <- Phihat %*% tcrossprod(invdel, Phihat)
      temp1 <- kronecker(m1, m2)
      d1 <- etahat %*% tcrossprod(invome, etahat)
      m3 <- kronecker(d1, Delta0hat)
      m4 <- kronecker(Deltahat, invdel0) + kronecker(invdel, 
           Delta0hat) - 2 * kronecker(diag(q), diag(p - q))
      m5 <- m3 + m4
      invm5 <- chol2inv(chol(m5))
      m6 <- kronecker(tcrossprod(Gammahat, etahat), Phi0hat)
      temp2 <- m6 %*% tcrossprod(invm5, m6)
      d2 <- crossprod(etahat, Deltahat) %*% etahat
      m8 <- kronecker(invome0, d2) + kronecker(invome0, Omegahat)
      m9 <- kronecker(Omega0hat, invome) - 2 * kronecker(diag(r - u), diag(u))
      m10 <- m8 + m9
      invm10 <- chol2inv(chol(m10))
      m11 <- kronecker(Gamma0hat, Phihat %*% etahat)
      temp3 <- m11 %*% tcrossprod(invm10, m11)
      covMatrix <- temp1 + temp2 + temp3
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- asyFm/asySE
    }
  }

  return(list(beta = betahat, SigmaYcX = SigmaYcXhat, SigmaX = SigmaXhat, 
              Gamma = Gammahat, Gamma0 = Gamma0hat,
              eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, 
              mu = muhat, Phi = Phihat, Phi0 = Phi0hat,
              Delta = Deltahat, Delta0 = Delta0hat,
              loglik = loglik, covMatrix = covMatrix, 
              asySE = asySE, ratio = ratio, n = n))
}
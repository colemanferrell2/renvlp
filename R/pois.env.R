pois.env <- function(X, Y, u, asy = TRUE, init = NULL) {
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
  
  sigY <- stats::cov(Y) * (n - 1)/n
  sigYX <- stats::cov(Y, X) * (n - 1)/n
  sigX <- stats::cov(X) * (n - 1)/n
  invsigY <- chol2inv(chol(sigY))
  tmp <- pois.envMU(X, Y, u)
  
  if (!is.null(init)) {
      if (nrow(init) != p || ncol(init) != u) stop("The initial value should have p rows and u columns.")
      tmp0 <- qr.Q(qr(init), complete = TRUE)
      tmp$Gammahat <- tmp0[, 1:u]
      tmp$Gamma0hat <- tmp0[, (u+1):p]
      GX <- X %*% tmp$Gammahat
      fit <- stats::glm(Y ~ GX, family = stats::poisson)
      tmp$muhat <- as.vector(fit$coefficients[1])
      tmp$etahat <- matrix(fit$coefficients[2 : (u + 1)])
      beta <- tmp$Gammahat %*% tmp$etahat
      theta <- matrix(1, n, 1) %*% tmp$muhat + X %*% beta
      c.theta <- - exp(theta)
      c.theta.mean <- mean(c.theta)
      tmp$weighthat <- c.theta / c.theta.mean
      tmp$Vhat <- theta + ((Y - c.theta) / tmp$weighthat)
      W <- diag(as.vector(tmp$weighthat))
      wx <- W %*% X
      mean.wx <- apply(wx, 2, mean)
      wxx <- X - matrix(1, nrow = n) %*% mean.wx
      sigwx <- crossprod(wxx, W) %*% wxx / n
      wv <- W %*% tmp$Vhat
      mean.wv <- mean(wv)
      wvv <- tmp$Vhat - matrix(1, nrow = n) %*% mean.wv
      sigwxv <- crossprod(wxx, W) %*% wvv / n
      inv.sigwx <- chol2inv(chol(sigwx))
      tmp$avar <- inv.sigwx / (- c.theta.mean)
      M <- sigX
      MU <- sigX
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2,
      1/tmp.MU$values, "*") %*% t(tmp.MU$vectors)
      e1 <- eigen(t(tmp$Gammahat) %*% M %*% tmp$Gammahat)
      e2 <- eigen(t(tmp$Gammahat) %*% invMU %*% tmp$Gammahat)
      e3 <- eigen(M)
      e4 <- matrix(1, n, 1) %*% tmp$muhat + X %*% tmp$Gammahat %*% tmp$etahat
      e5 <- crossprod(Y, e4) - colSums(exp(e4))
      tmp$objfun <- - n/2 * (sum(log(e1$values)) + sum(log(e2$values)) + sum(log(e3$values))) + e5
  }
  
  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  objfun <- tmp$objfun
  muhat <- tmp$muhat
  etahat <- tmp$etahat
  weighthat <- tmp$weighthat
  Vhat <- tmp$Vhat
  avarhat <- tmp$avar
  
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL
  if (u == 0) {
    Omegahat <- NULL
    Omega0hat <- sigX
    betahat <- matrix(0, p, 1)
    SigmaXhat <- sigX
    loglik <- objfun
    if (asy == T) 
      ratio <- matrix(1, p, r)
  }
  else if (u == p) {
    Omegahat <- sigX
    Omega0hat <- NULL
    betahat <- etahat
    SigmaXhat <- sigX
    loglik <- objfun
    if (asy == T) {
      covMatrix <- avarhat
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- matrix(1, p, r)
    }
  }
  else {
    invsigX <- chol2inv(chol(sigX))
    Omegahat <- crossprod(Gammahat, sigX) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigX) %*% Gamma0hat
    betahat <- Gammahat %*% etahat
    SigmaXhat <- Gammahat %*% tcrossprod(Omegahat, 
                 Gammahat) + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
    PGamma <- tcrossprod(Gammahat)
    loglik <- objfun
    if (asy == T) {
      options(warn = -1) 
      fit <- stats::glm(Y ~ X, family = stats::poisson)
      mu <- as.vector(fit$coefficients[1])
      beta <- as.vector(fit$coefficients[2 : (p + 1)])
      theta <- matrix(1, n, 1) %*% mu + X %*% beta
      c.theta <- - exp(theta)
      c.theta.mean <- mean(c.theta)
      weight <- c.theta / c.theta.mean
      V <- theta + ((Y - c.theta) / weight)
      W <- diag(as.vector(weight))
      wx <- W %*% X
      mean.wx <- apply(wx, 2, mean)
      wxx <- X - matrix(1, nrow = n) %*% mean.wx
      sigwx <- crossprod(wxx, W) %*% wxx / n
      wv <- W %*% V
      mean.wv <- mean(wv)
      wvv <- V - matrix(1, nrow = n) %*% mean.wv
      sigwxv <- crossprod(wxx, W) %*% wvv / n
      inv.sigwx <- chol2inv(chol(sigwx))
      avar <- inv.sigwx / (- c.theta.mean)
      covMatrix <- avar
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = p)
      PGamma <- tcrossprod(Gammahat)
      t1 <- avarhat
      T1 <- PGamma %*% t1 %*% PGamma
      invt1 <- chol2inv(chol(t1))
      t2 <- kronecker(etahat, t(Gamma0hat))
      invOme <- chol2inv(chol(Omegahat))
      invOme0 <- chol2inv(chol(Omega0hat))
      t3 <- kronecker(Omegahat, invOme0) + kronecker(invOme, 
            Omega0hat) - 2 * kronecker(diag(u), diag(p - u))
      M <- t2 %*% tcrossprod(invt1, t2) + t3
      invM <- chol2inv(chol(M))
      T2 <- crossprod(t2, invM) %*% t2
      covMatrix <- T1 + T2
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = p)
      ratio <- asyFm/asySE
    }
  }
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat,
              beta = as.matrix(betahat), SigmaX = SigmaXhat, eta = etahat, 
              Omega = Omegahat, Omega0 = Omega0hat, 
              loglik = loglik, n = n, covMatrix = covMatrix, asySE = asySE, 
              ratio = ratio))
}
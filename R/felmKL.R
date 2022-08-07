felmKL <- function(X, Y, ux, uy, t1, t2, knots = c(0, 0.25, 0.5, 0.75, 1)) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)  
  dimX <- dim(X)
  n <- dimX[1]
  
  if (min(knots) < 0) stop("Location of the knots cannot be negative.")
  if (is.null(knots)) {
    knots <- orthogonalsplinebasis::expand.knots(c(0, 0.25, 0.5, 0.75, 1), order = 4)
  } else {
    knots <- knots / max(knots)
    knots <- orthogonalsplinebasis::expand.knots(knots, order = 4)
  }
  
  if (nrow(Y)!= nrow(X)) stop("X and Y should have the same number of observations.")
  if (min(t1) < 0) stop("Time points cannot be negative.")
  if (min(t2) < 0) stop("Time points cannot be negative.")
  if (max(t1) > 1) t1 <- t1 / max(t1)
  if (max(t2) > 1) t2 <- t2 / max(t2)
  
  spbasis <- orthogonalsplinebasis::SplineBasis(knots, order = 4)
  Gx <- orthogonalsplinebasis::GramMatrix(spbasis)
  eigX <- eigen(Gx)
  Gx12 <- eigX$vectors %*% diag(sqrt(eigX$values)) %*% t(eigX$vectors)
  dx <- dim(Gx)[1]
  
  X.cord <- matrix(0, n, dx)
  Y.cord <- matrix(0, n, dx)
  basis.value.X <- orthogonalsplinebasis::evaluate(spbasis, t1)
  basis.value.Y <- orthogonalsplinebasis::evaluate(spbasis, t2)
  
  for (i in 1:n) {
      X.cord[i,] <- stats::lm(X[i, ] ~ basis.value.X -1)$coefficients
      Y.cord[i,] <- stats::lm(Y[i, ] ~ basis.value.Y -1)$coefficients
  }

  tmpX <- eigen(Gx12 %*% stats::cov(X.cord) %*% t(Gx12))
  tmpY <- eigen(Gx12 %*% stats::cov(Y.cord) %*% t(Gx12))
  X.melm <- X.cord %*% Gx12 %*% tmpX$vectors
  Y.melm <- Y.cord %*% Gx12 %*% tmpY$vectors
  Gx12inv <- solve(Gx12)
  phihat.cord <- Gx12inv %*% tmpX$vectors
  psihat.cord <- Gx12inv %*% tmpY$vectors
  
  m <- stenv(X.melm, Y.melm, ux, uy)
  mfull <- stenv(X.melm, Y.melm, dx, dx)
  
  psihat.Y <- basis.value.Y %*% psihat.cord
  
  fitted.env <- (X.melm %*% m$beta + matrix(1, n, 1) %*% t(m$mu)) %*% t(psihat.Y)
  fitted.full <- (X.melm %*% mfull$beta + matrix(1, n, 1) %*% t(mfull$mu)) %*% t(psihat.Y)
  
  return(list(beta = m$beta, betafull = mfull$beta, alpha = m$mu, alphafull = mfull$mu, phihat.cord = phihat.cord, psihat.cord = psihat.cord, fitted.env = fitted.env, fitted.full = fitted.full))
}

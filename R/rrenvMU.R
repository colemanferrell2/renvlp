rrenvMU <- function(M, U, u, d) {

  dimM <- dim(M)
  dimU <- dim(U)
  r <- dimM[1]
  
  if (dimM[1] != dimM[2] & dimU[1] != dimU[2]) stop("M and U should be square matrices.")
  if (dimM[1] != dimU[1]) stop("M and U should have the same dimension.")
  if (qr(M, tol = 1e-10)$rank < r) stop("M should be positive definite.")
  if (u > r & u <= d) stop("u should be between d and r.")
  

	if (u == 0) {
	
		Gammahat <- NULL
		Gamma0hat <- diag(r)
		MU <- M + U
		tmp.MU <- eigen(MU)
		objfun <- sum(log(tmp.MU$values))
		
	} else if (u == r) {
	  
	  Gammahat <- diag(r)
	  Gamma0hat <- NULL
	  tmp.M <- eigen(M)
	  objfun <- sum(log(tmp.M$values))
	  
	} else if (u == r - 1) {
	  
	  maxiter = 100
	  ftol = 1e-3
	  
	  MU <- M + U
	  tmp.MU <- eigen(MU)
	  invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
	  invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
	  
	  midmatrix <- U
	  startv <- function(a) t(a) %*% midmatrix %*% a
	  tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
	  tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
	  init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
	  
#	  if (qr(MU)$rank == r) {
	    eig1 <- eigen(t(init) %*% M %*% init)
	    eig2 <- eigen(t(init) %*% invMU %*% init)
		  GMGnhalf <- sweep(eig1$vectors, MARGIN = 2, sqrt(1 / eig1$values), '*') %*% t(eig1$vectors)	
		  a <- eigen(MU)
		  cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
		  obj1 <- sum(log(eig1$values)) + sum(log(eig2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))
	    
	    midmatrix <- invMU2 %*% tcrossprod(U, invMU2) 
	    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
	    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
	    init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
	    e1 <- eigen(t(init.MU) %*% M %*% init.MU)
	    e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
		  GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)
		  a <- eigen(MU)
		  cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
		  obj2 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))		
	    if (obj2 < obj1) {
	      init <- init.MU
	      obj1 <- obj2
	    }
	    
#	    if (qr(M)$rank == r) {
	      tmp.M <- eigen(M)
	      midmatrix <- U
	      tmp2.M <- apply(tmp.M$vectors, 2, startv)
	      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
	      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
	      e1 <- eigen(t(init.M) %*% M %*% init.M)
	      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
		    GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)
		    a <- eigen(MU)
		    cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
		    obj3 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))			
	      if (obj3 < obj1) {
	        init <- init.M
	        obj1 <- obj3
	      }
	      
	      invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
	      midmatrix <- invM2 %*% tcrossprod(U, invM2) 
	      tmp2.M <- apply(tmp.M$vectors, 2, startv)
	      tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
	      init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				
	      
	      e1 <- eigen(t(init.M) %*% M %*% init.M)
	      e2 <- eigen(t(init.M) %*% invMU %*% init.M)
		    GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)			
		    a <- eigen(MU)
		    cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
		    obj4 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))			
	      if (obj4 < obj1) {
	        init <- init.M
	        obj1 <- obj4
	      }

	  
	  GEidx <- GE(init)
	  Ginit <- init %*% solve(init[GEidx[1:u], ])		
	  
	  j <- GEidx[r]
	  
	  g <- as.matrix(Ginit[j, ])
	  t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
	  t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
	  t6 <- crossprod(Ginit[-j, ], as.matrix(MU[-j, j])) / MU[j, j]
	  
	  GUGt2 <- g + t2
	  GUG <- crossprod(Ginit, (M %*% Ginit)) - tcrossprod(GUGt2, GUGt2) * M[j, j]
	  
	  GVGt2 <- g + t3
	  GVG <- crossprod(Ginit, (invMU %*% Ginit)) - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
	  
	  GWGt2 <- g + t6
	  GWG <- crossprod(Ginit, (MU %*% Ginit)) - tcrossprod(GWGt2, GWGt2) * MU[j, j]
      GWG <- as.matrix(Matrix::nearPD(GWG)$mat)
	  
	  invC1 <- chol2inv(chol(GUG))
	  invC2 <- chol2inv(chol(GVG))
	  invC6 <- chol2inv(chol(GWG))
	  
	  
	  fobj <- function(x) {
	    tmp2 <- x + t2
	    tmp3 <- x + t3
	    tmp6 <- x + t6
	    T2 <- invC1 %*% tmp2	
	    T3 <- invC2 %*% tmp3
	    T6 <- invC6 %*% tmp6
	    H1 <- tcrossprod(tmp2) * M[j, j] + GUG
	    H1 <- as.matrix(Matrix::nearPD(H1)$mat)
	    H2 <- tcrossprod(tmp6) * MU[j, j] + GWG
	    a <- eigen(H1)
	    tmp4 <- a$vectors %*% diag(sqrt(1 / a$values)) %*% t(a$vectors)
	    tmp5 <- crossprod(tmp4, H2) %*% tmp4
	    tmp5 <- as.matrix(Matrix::nearPD(tmp5)$mat)
	    b <- eigen(tmp5)
	    y <- -2 * log(1 + sum(x^2)) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3)) + sum(log(b$values[(d+1):u]))	
	    return(y)
	  }
	  
	  gobj <- function(x) {
	    tmp2 <- x + t2
	    tmp3 <- x + t3
	    tmp6 <- x + t6
	    T2 <- invC1 %*% tmp2	
	    T3 <- invC2 %*% tmp3
	    T6 <- invC6 %*% tmp6
	    H1 <- tcrossprod(tmp2) * M[j, j] + GUG
	    H2 <- tcrossprod(tmp6) * MU[j, j] + GWG
	    a <- eigen(H1)
	    tmp4 <- sweep(a$vectors, MARGIN = 2, sqrt(1 / a$values), '*') %*% t(a$vectors)
	    tmp5 <- crossprod(tmp4, H2) %*% tmp4
	    b <- eigen(tmp5)
	    H1inv <- sweep(a$vectors, MARGIN = 2, 1 / a$values, '*') %*% t(a$vectors)
	    U <- tmp4 %*% b$vectors
	    V <- sweep(a$vectors, MARGIN = 2, sqrt(a$values), '*') %*% t(a$vectors) %*% b$vectors
	    UVtmp <- sweep(as.matrix(V[, (d+1):u]), MARGIN = 2, 1 / b$values[(d+1):u], '*') %*% t(U[, (d+1):u])
	    dlamdZ <- matrix(UVtmp, nrow = 1)
	    dh1da <- M[j, j] * (kronecker(tmp2, diag(u)) + kronecker(diag(u), tmp2))
	    dh2da <- invMU[j, j] * (kronecker(tmp3, diag(u)) + kronecker(diag(u), tmp3))
	    dlamda <- dlamdZ %*% (-kronecker(H2 %*% H1inv, H1inv) %*% dh1da + kronecker(diag(u), H1inv) %*% dh2da) 
	    y <- -4 	* x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))	+ t(dlamda)	
	    return(y)
	  }
	  
	  i <- 1
	  while (i < maxiter) {
	    
		if (sum(Re(gobj(Ginit[j,]))^2) < 1e10) {
            res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
			if (sum(res$par^2) < 1e10) Ginit[j, ] <- res$par
		}
	    
	    a <- qr(Ginit)
	    Gammahat <- qr.Q(a)
	    e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
	    e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)				
		GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)			
		gmg <- GMGnhalf %*% t(Gammahat) %*% MU %*% Gammahat %*% GMGnhalf
		gmg <- as.matrix(Matrix::nearPD(gmg)$mat)
		obj5 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(gmg)$values[(d+1):u]))	
	    if (abs(obj1 - obj5) < ftol * abs(obj1)) {
	      break
	    } else {
	      obj1 <- obj5
	      i <- i + 1
	    }
	  }
	  
	  Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r, drop = FALSE]
	  objfun <- obj5 + sum(log(tmp.MU$values))
		
	} else {
    
		maxiter = 100
		ftol = 1e-3
		

		MU <- M + U
		tmp.MU <- eigen(MU)
		invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
    invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
		
    midmatrix <- U
    startv <- function(a) t(a) %*% midmatrix %*% a
		tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
		tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
		init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
			
#		if (qr(MU)$rank == r) {
			eig1 <- eigen(t(init) %*% M %*% init)
			eig2 <- eigen(t(init) %*% invMU %*% init)
			GMGnhalf <- sweep(eig1$vectors, MARGIN = 2, sqrt(1 / eig1$values), '*') %*% t(eig1$vectors)
			a <- eigen(MU)
			cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
			obj1 <- sum(log(eig1$values)) + sum(log(eig2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))
			
			midmatrix <- invMU2 %*% tcrossprod(U, invMU2) 
			tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
			tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
			init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
			e1 <- eigen(t(init.MU) %*% M %*% init.MU)
			e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
			GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)
			a <- eigen(MU)
			cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
			obj2 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))		
			if (obj2 < obj1) {
				init <- init.MU
				obj1 <- obj2
			}
			
#			if (qr(M)$rank == r) {
				tmp.M <- eigen(M)
				midmatrix <- U
				tmp2.M <- apply(tmp.M$vectors, 2, startv)
				tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
				init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
				e1 <- eigen(t(init.M) %*% M %*% init.M)
				e2 <- eigen(t(init.M) %*% invMU %*% init.M)
				GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)
				a <- eigen(MU)
				cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
				obj3 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))				
				if (obj3 < obj1) {
					init <- init.M
					obj1 <- obj3
				}
				
				invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
				midmatrix <- invM2 %*% tcrossprod(U, invM2) 
				tmp2.M <- apply(tmp.M$vectors, 2, startv)
				tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
				init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				

				e1 <- eigen(t(init.M) %*% M %*% init.M)
				e2 <- eigen(t(init.M) %*% invMU %*% init.M)
				GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)			
				a <- eigen(MU)
				cc <- GMGnhalf %*% t(init) %*% a$vector %*% diag(sqrt(a$values))
				obj4 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(tcrossprod(cc))$values[(d+1):u]))			
				if (obj4 < obj1) {
					init <- init.M
					obj1 <- obj4
				}
#			}
#		}
		
		GEidx <- GE(init)
		Ginit <- init %*% solve(init[GEidx[1:u], ])

		GUG <- crossprod(Ginit, (M %*% Ginit))	
		GVG <- crossprod(Ginit, (invMU %*% Ginit))
		GWG <- crossprod(Ginit, (MU %*% Ginit))		

		
		t4 <- crossprod(Ginit[GEidx[(u+1):r],], Ginit[GEidx[(u+1):r], ]) + diag(u)
		i <- 1
		while (i < maxiter) {
			for (j in GEidx[(u+1):r]) {
				g <- as.matrix(Ginit[j, ])
				t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
				t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
				t6 <- crossprod(Ginit[-j, ], as.matrix(MU[-j, j])) / MU[j, j]
				
				GUGt2 <- g + t2
				GUG <- GUG - tcrossprod(GUGt2, GUGt2) * M[j, j]
				GUG <- as.matrix(Matrix::nearPD(GUG)$mat)
				
				GVGt2 <- g + t3
				GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
				
				GWGt2 <- g + t6
				GWG <- GWG - tcrossprod(GWGt2, GWGt2) * MU[j, j] 
				GWG <- as.matrix(Matrix::nearPD(GWG)$mat)
				
				t4 <- t4 - tcrossprod(g, g)
				invC1 <- chol2inv(chol(GUG))
				invC2 <- chol2inv(chol(GVG))
				invt4 <- chol2inv(chol(t4))		
				invC6 <- chol2inv(chol(GWG))			
									
									
				fobj <- function(x) {
					tmp2 <- x + t2
					tmp3 <- x + t3
					tmp6 <- x + t6
					T1 <- invt4 %*% x
					T2 <- invC1 %*% tmp2	
					T3 <- invC2 %*% tmp3
					T6 <- invC6 %*% tmp6
					H1 <- tcrossprod(tmp2) * M[j, j] + GUG
					H1 <- as.matrix(Matrix::nearPD(H1)$mat)
					H2 <- tcrossprod(tmp6) * MU[j, j] + GWG
					a <- eigen(H1)
					tmp4 <- sweep(a$vectors, MARGIN = 2, sqrt(1 / a$values), '*') %*% t(a$vectors)
					tmp5 <- crossprod(tmp4, H2) %*% tmp4
					b <- eigen(tmp5)
					y <- -2 * log(1 + x %*% T1) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3)) + sum(log(b$values[(d+1):u]))
					return(y)
				}
			
				gobj <- function(x) {
					tmp2 <- x + t2
					tmp3 <- x + t3
					tmp6 <- x + t6
					T1 <- invt4 %*% x
					T2 <- invC1 %*% tmp2	
					T3 <- invC2 %*% tmp3
					T6 <- invC6 %*% tmp6
					H1 <- tcrossprod(tmp2) * M[j, j] + GUG
					H2 <- tcrossprod(tmp6) * MU[j, j] + GWG
					a <- eigen(H1)
					tmp4 <- sweep(a$vectors, MARGIN = 2, sqrt(1 / a$values), '*') %*% t(a$vectors)
					tmp5 <- crossprod(tmp4, H2) %*% tmp4
					b <- eigen(tmp5)
					H1inv <- sweep(a$vectors, MARGIN = 2, 1 / a$values, '*') %*% t(a$vectors)
					U <- tmp4 %*% b$vectors
					V <- sweep(a$vectors, MARGIN = 2, sqrt(a$values), '*') %*% t(a$vectors) %*% b$vectors
					UVtmp <- sweep(as.matrix(V[, (d+1):u]), MARGIN = 2, 1 / b$values[(d+1):u], '*') %*% t(U[, (d+1):u])
					dlamdZ <- matrix(UVtmp, nrow = 1)
					dh1da <- M[j, j] * (kronecker(tmp2, diag(u)) + kronecker(diag(u), tmp2))
					dh2da <- invMU[j, j] * (kronecker(tmp3, diag(u)) + kronecker(diag(u), tmp3))
					dlamda <- dlamdZ %*% (-kronecker(H2 %*% H1inv, H1inv) %*% dh1da + kronecker(diag(u), H1inv) %*% dh2da) 
					y <- -4 	* T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))	+ t(dlamda)
					return(y)
				}

				if (sum(Re(gobj(Ginit[j,]))^2) < 1e10) {
                    res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
					if (sum(res$par^2) < 1e10) Ginit[j, ] <- res$par
				}
				g <- as.matrix(Ginit[j, ])
				t4 <- t4 + tcrossprod(g, g)
				GUGt2 <- g + t2
				GUG <- GUG + tcrossprod(GUGt2, GUGt2) * M[j, j]
				
				GVGt2 <- g + t3
				GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
				
				GWGt2 <- g + t6
				GWG <- GWG + tcrossprod(GWGt2, GWGt2) * MU[j, j]
				
				
			}
			a <- qr(Ginit)
			Gammahat <- qr.Q(a)
			e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
			e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)		
			GMGnhalf <- sweep(e1$vectors, MARGIN = 2, sqrt(1 / e1$values), '*') %*% t(e1$vectors)			
			obj5 <- sum(log(e1$values)) + sum(log(e2$values)) + sum(log(eigen(GMGnhalf %*% t(Gammahat) %*% MU %*% Gammahat %*% GMGnhalf)$values[(d+1):u]))	
			if (abs(obj1 - obj5) < ftol * abs(obj1)) {
				break
			} else {
				obj1 <- obj5
				i <- i + 1
			}
		}

		Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r]
		objfun <- obj5 + sum(log(tmp.MU$values))
		
	}
	return(list(Gammahat = Gammahat, Gamma0hat = Gamma0hat, objfun = objfun))
}
	




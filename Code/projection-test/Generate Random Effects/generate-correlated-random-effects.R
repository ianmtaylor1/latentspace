library(MASS)

n <- 20
reps.per.X <- 5
num.X <- 5
alpha.list <- seq(0,1,by=.1)
sigma.list <- 1:5

########################################

# Raise a matrix to an integer power (inefficient, but okay)
matpow <- function(M, p) {
  res <- M
  for (i in 2:p) {
    res <- res %*% M
  }
  return(res)
}

# Do the eta-prime transformation
eta.prime <- function(alpha, sigma, eta, I, P) {
  # Divide eta into components
  eta.orth <- (I - P) %*% eta
  eta.par <- P %*% eta
  # Reassemble and return
  (sqrt(1-alpha) * eta.orth / sqrt(sum(eta.orth^2)) + 
      sqrt(alpha) * eta.par / sqrt(sum(eta.par^2))) * sigma
}

########################################

results <- NULL

for (gen in 1:num.X) {
  X <- cbind(1, matrix(rnorm(4*n^2), ncol=4, nrow=n^2))
  
  A <- matrix(0, ncol=n, nrow=n^2)
  B <- matrix(0, ncol=n, nrow=n^2)
  for (i in 1:n) {
    Y <- matrix(0, nrow=n, ncol=n)
    Y[i,] <- 1
    A[,i] <- c(Y)
    B[,i] <- c(t(Y))
  }
  AB <- cbind(A,B)
  
  # Various projection matrices onto ...
  PX <- X %*% solve(t(X) %*% X) %*% t(X)     # The column space of X
  PA <- A %*% solve(t(A) %*% A) %*% t(A)     # The column space of A
  PB <- B %*% solve(t(B) %*% B) %*% t(B)     # The column space of B
  PAB <- AB %*% ginv(t(AB) %*% AB) %*% t(AB) # The column space of AB 
  I <- diag(n^2)                             # R^(n^2)
  PAX <- matpow(PA %*% PX, 100)              # (C(A) intersect C(X))
  PBX <- matpow(PB %*% PX, 100)              # (C(B) intersect C(X))
  PABX <- matpow(PAB %*% PX, 100)            # (C(AB) intersect C(X))
  PXA <- (PA %*% X) %*% ginv(t(X) %*% PA %*% X) %*% t(PA %*% X) # The column space of PA%*%X
  PXB <- (PB %*% X) %*% ginv(t(X) %*% PB %*% X) %*% t(PB %*% X) # The column space of PB%*%X
  PXAB <- (PAB %*% X) %*% ginv(t(X) %*% PAB %*% X) %*% t(PAB %*% X) # The column space of PAB%*%X
  
  for (rep in 1:reps.per.X) {
    eta <- A %*% rnorm(n)
    
    for (alpha in alpha.list) {
      for (sigma in sigma.list) {
        eta.0 <- eta.prime(alpha, sigma, eta, I, PX)
        eta.1 <- PA %*% eta.0
        eta.2 <- eta.prime(alpha, sigma, eta, PAB, PABX)
        eta.3 <- eta.prime(alpha, sigma, eta, PAB, PXAB)
        
        results <- rbind(results,
                         list(alpha=alpha, sigma=sigma, rep=rep, Xgen=gen,
                              propvar.0=sum((PX%*%eta.0)^2)/sum(eta.0^2),
                              propvar.1=sum((PX%*%eta.1)^2)/sum(eta.1^2),
                              propvar.2=sum((PX%*%eta.2)^2)/sum(eta.2^2),
                              propvar.3=sum((PX%*%eta.3)^2)/sum(eta.3^2),
                              mag.0=sqrt(sum(eta.0^2)),
                              mag.1=sqrt(sum(eta.1^2)),
                              mag.2=sqrt(sum(eta.2^2)),
                              mag.3=sqrt(sum(eta.3^2)))
        )
      }
    }
  }
}





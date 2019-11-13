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

genX <- function(n) {
  #### Generate design matrix X for dyad covariates
  # Start empty
  Xd <- array(0, dim=c(n,n,4) )
  # Fill first covariate with binary indicator of shared class
  groups <- sample(0:1, n, TRUE) * 2 - 1
  if (all(groups == 1) || all(groups == -1)) {
    k <- sample(1:n, 1)[1]
    groups[k] <- -groups[k]
  }
  Xd[,,1] <- outer(groups, groups, "*")
  # Fill second covariate with "distances" between normally distributed coordinates
  locations <- rnorm(n=n, mean=0, sd=1)
  dist <- function(x,y) {
    abs(x - y)
  }
  Xd[,,2] <- outer(locations, locations, "dist")
  # Fill third covariate with all random standard normal draws
  Xd[,,3] <- matrix(rnorm(n=n*n, mean=0, sd=1), nrow=n, ncol=n)
  # Fill fourth covariate with a random row effect
  Xd[,,4] <- matrix(rnorm(n=n, mean=0, sd=1), nrow=n, ncol=n)
  # Set all diagonals to NA
  #for (i in 1:n) {
  #  for (k in 1:length(beta)) {
  #    Xd[i,i,k] <- NA
  #  }
  #}
  
  dimnames(Xd)[[3]] <- c("X1","X2","X3","X4")
  
  mX <- apply(Xd,3,c)
  
  return(cbind(1, mX))
}

########################################

results <- NULL

for (gen in 1:num.X) {
  X <- genX(n)
  
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
  I <- diag(n^2)                             #R^(n^2)
  PX <- X %*% solve(t(X) %*% X) %*% t(X)     # The column space of X
  PAB <- AB %*% ginv(t(AB) %*% AB) %*% t(AB) # The column space of AB 
  PABX <- matpow(PAB %*% PX, 200)            # (C(AB) intersect C(X))
  
  for (rep in 1:reps.per.X) {
    eta <- AB %*% rnorm(2*n)
    
    for (alpha in alpha.list) {
      for (sigma in sigma.list) {
        eta.0 <- eta.prime(alpha, sigma, eta, I, PX)
        eta.2 <- eta.prime(alpha, sigma, eta, PAB, PABX)
        
        results <- rbind(results,
                         list(alpha=alpha, sigma=sigma, rep=rep, Xgen=gen,
                              propvar.0=sum((PX%*%eta.0)^2)/sum(eta.0^2),
                              propvar.2=sum((PX%*%eta.2)^2)/sum(eta.2^2),
                              mag.0=sqrt(sum(eta.0^2)),
                              mag.2=sqrt(sum(eta.2^2)),
                              rhoX1.0=cor(X[,2], eta.0),
                              rhoX2.0=cor(X[,3], eta.0),
                              rhoX3.0=cor(X[,4], eta.0),
                              rhoX4.0=cor(X[,5], eta.0),
                              rhoX1.2=cor(X[,2], eta.2),
                              rhoX2.2=cor(X[,3], eta.2),
                              rhoX3.2=cor(X[,4], eta.2),
                              rhoX4.2=cor(X[,5], eta.2))
        )
      }
    }
  }
}





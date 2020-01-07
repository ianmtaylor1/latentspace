library(MASS)
library(foreach)

outfile <- "Gen-rnd-effect-correlation-plots.pdf"

n <- 20
reps.per.X <- 3
num.X <- 4
alpha.list <- seq(0,1,by=.1)
sigma.list <- 1:2

########################################

# Do the eta-prime transformation
eta.prime <- function(alpha, sigma, eta, P, P.orth) {
  n <- length(eta)
  # Divide eta into components
  eta.orth <- P.orth %*% eta
  eta.par <- P %*% eta
  # Reassemble and return
  (sqrt(1-alpha) * eta.orth / sqrt(sum(eta.orth^2)) + 
      sqrt(alpha) * eta.par / sqrt(sum(eta.par^2))) * sigma * n
}

genX <- function(n) {
  #### Generate design matrix X for dyad covariates
  # Start empty
  Xd <- array(0, dim=c(n,n,3) )
  # Fill first covariate with a random row effect
  Xd[,,1] <- matrix(rnorm(n=n, mean=0, sd=1), nrow=n, ncol=n)
  # Fill the second covariate with an "almost" row covariate
  Xd[,,2] <- matrix(rnorm(n=n, mean=0, sd=1), nrow=n, ncol=n)
  Xd[,,2] <- Xd[,,2] + matrix(rnorm(n=n^2, mean=0, sd=1/n), nrow=n, ncol=n)
  # Fill third covariate with all random standard normal draws
  Xd[,,3] <- matrix(rnorm(n=n*n, mean=0, sd=1), nrow=n, ncol=n)
  
  dimnames(Xd)[[3]] <- c("X1","X2","X3")
  
  mX <- apply(Xd,3,c)
  
  return(cbind(1, mX))
}

########################################

results <- foreach(gen = 1:num.X, .combine="rbind") %do% {
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
  # C(PAB %*% X) - the column space of X projected onto AB.
  PABX <- (PAB %*% X) %*% ginv(t(X) %*% PAB %*% X) %*% (t(X) %*% PAB)
  # Attempting to do a "weighted" projection
  PABXw <- (PAB %*% X) %*% diag(1 / diag(t(X) %*% X)) %*% (t(X) %*% PAB)
  
  foreach(rep = 1:reps.per.X, .combine="rbind") %do% {
    eta <- AB %*% rnorm(2*n)
    
    foreach(alpha = alpha.list, .combine="rbind") %do% {
      foreach(sigma = sigma.list, .combine="rbind") %do% {
        eta.4 <- eta.prime(alpha, sigma, eta, PABXw, PAB - PABX)
        
        data.frame(list(alpha=alpha, sigma=sigma, rep=rep, Xgen=gen,
                   propvar.4=sum((PX%*%eta.4)^2)/sum(eta.4^2),
                   cancor.4=cancor(X, eta.4)$cor,
                   mag.4=sqrt(sum(eta.4^2)),
                   rhoX1.4=cor(X[,2], eta.4),
                   rhoX2.4=cor(X[,3], eta.4),
                   rhoX3.4=cor(X[,4], eta.4)))
        
      }
    }
  }
}
results <- as.data.frame(results)

pdf(outfile, width=8, height=6)

par(mfrow=c(2,2))
# Canonical correlation
for (i in 1:num.X) {
  plot(x=results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'alpha'], 
       y=results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'cancor.4'],
       main=paste("Canonical correlation (X #",i,")",sep=""), ylim=c(0,1),
       xlab=expression(alpha), ylab="Canonical correlation")
}
# (Absolute) Correlation with the row covariate
for (i in 1:num.X) {
  plot(x=results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'alpha'], 
       y=abs(results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'rhoX1.4']),
       main=paste("Correlation with row covariate (X #",i,")",sep=""), ylim=c(0,1),
       xlab=expression(alpha), ylab=expression(abs(rho[1])))
}
# (Absolute) Correlation with the almost row covariate
for (i in 1:num.X) {
  plot(x=results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'alpha'], 
       y=abs(results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'rhoX2.4']),
       main=paste("Correlation with almost row covariate (X #",i,")",sep=""), ylim=c(0,1),
       xlab=expression(alpha), ylab=expression(abs(rho[2])))
}
# (Absolute) Correlation with dyad covariate
for (i in 1:num.X) {
  plot(x=results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'alpha'], 
       y=abs(results[(results$sigma==1)&(results$Xgen==i)&(results$rep==1), 'rhoX3.4']),
       main=paste("Correlation with unstructured covariate (X #",i,")",sep=""), ylim=c(0,1),
       xlab=expression(alpha), ylab=expression(abs(rho[3])))
}


dev.off()


library(amenhs)
library(foreach)
library(MASS)

n <- 20       # Number of nodes in the network

runs <- 5
reps <- 10

mcmc.iter <- 10000 # Number of iterations of the MCMC algorithm to do (after "burn") (default: 10000)

intercept <- -1
beta <- c(1,1,1) # True coefficients (not including intercept)
gamma <- 1       # Coefficient for unobserved data

# ame() resets random seed, so seed each loop with new number
seedbase <- 20200107

######################################################

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
  
  return(Xd)
  
  #mX <- apply(Xd,3,c)
  
  #return(cbind(1, mX))
}

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
PAB <- AB %*% ginv(t(AB) %*% AB) %*% t(AB) # The column space of AB 


#######################################

# For each run
outdata <- foreach(run=1:runs, .combine="rbind") %do% {
  cat("Run:",run,"\n")
  set.seed(seedbase + 10000*run)
  #### Generate design matrix X for dyad covariates
  Xd <- genX(n)
  mX <- cbind(1,apply(Xd,3,c))
  PX <- mX %*% solve(t(mX) %*% mX) %*% t(mX)     # The column space of X
  # C(PAB %*% X) - the column space of X projected onto AB.
  PABX <- (PAB %*% mX) %*% ginv(t(mX) %*% PAB %*% mX) %*% (t(mX) %*% PAB)
  # Attempting to do a "weighted" projection
  PABXw <- (PAB %*% mX) %*% diag(1 / diag(t(mX) %*% mX)) %*% (t(mX) %*% PAB)
  
  #### Generate unobserved covariate
  # ADDITIVE ROW AND COLUMN EFFECTS
  true.a <- rnorm(n=n)
  true.b <- rnorm(n=n)
  Z.additive <- outer(true.a, true.b, "+")
  
  Zc.low <-  eta.prime(0.0, 1, c(Z.additive), PABXw, PAB - PABX)
  Zc.high <- eta.prime(0.9, 1, c(Z.additive), PABXw, PAB - PABX)
  Z.low <-  matrix(Zc.low, nrow=n, ncol=n)
  Z.high <- matrix(Zc.high, nrow=n, ncol=n)
  
  diag(Z.additive) <- NA
  diag(Z.low) <- NA
  diag(Z.high) <- NA
  
  # For each rep
  foreach(rep=1:reps, .combine="rbind") %do% {
    cat("Rep:",rep,"-",strftime(Sys.time(),"%H:%M:%S"),"\n")
    set.seed(seedbase + 10000*run + rep)
    #### Simulate binary data with Probit
    eta <- intercept
    for (k in 1:length(beta)) {
      eta <- eta + beta[k]*Xd[,,k]
    }
    eta <- eta + matrix(rnorm(n = n*n), nrow=n, ncol=n)
    Y_binary_Zlow <- 1 * (eta + Z.low * gamma > 0) # Data generated with orthogonal unobserved covariate
    Y_binary_Zhigh <- 1 * (eta + Z.high * gamma > 0) # Data generated with correlated unobserved covariate
    
    # Fit a "projected" model for each Y
    foreach(fit=list(list(data=Y_binary_Zlow, alpha=0.0),
                     list(data=Y_binary_Zhigh, alpha=0.9)), 
            .combine="rbind") %do% {
              results <- amenhs::ame(Y=fit$data, Xdyad=Xd, family="bin",
                                     R=0, rvar=TRUE, cvar=TRUE, intercept=TRUE, 
                                     project=TRUE,
                                     plot=FALSE, gof=FALSE, print=FALSE)
              # Calculate useful things from the model
              foreach(col=c("intercept","X1.dyad","X2.dyad","X3.dyad"),
                      .combine="rbind") %do% {
                        list(run=run,
                             rep=rep,
                             alpha=fit$alpha,
                             col=col,
                             # Beta
                             beta.mean=mean(results$BETA[,col]),
                             beta.sd=sd(results$BETA[,col]),
                             beta.ci.low=quantile(results$BETA[,col],.025)[[1]],
                             beta.ci.high=quantile(results$BETA[,col],.975)[[1]],
                             # Delta
                             delta.mean=mean(results$DELTA[,col]),
                             delta.sd=sd(results$DELTA[,col]),
                             delta.ci.low=quantile(results$DELTA[,col],.025)[[1]],
                             delta.ci.high=quantile(results$DELTA[,col],.975)[[1]],
                             # Beta-tilde
                             betaT.mean=mean(results$BETA.TILDE[,col]),
                             betaT.sd=sd(results$BETA.TILDE[,col]),
                             betaT.ci.low=quantile(results$BETA.TILDE[,col],.025)[[1]],
                             betaT.ci.high=quantile(results$BETA.TILDE[,col],.975)[[1]],
                             # Mean square differences
                             msd.delta=mean((results$BETA[,col] - results$DELTA[,col])^2),
                             msd.tilde=mean((results$BETA[,col] - results$BETA.TILDE[,col])^2),
                             msd.dt=mean((results$DELTA[,col] - results$BETA.TILDE[,col])^2))
      }
    }
  }
}

write.csv(outdata, "2020-01-07 Quick Projection Results - Correlation.csv", row.names=FALSE)


###############################################################################
########## JUST DO ONE TO PLOT STUFF ##########################################
###############################################################################

set.seed(seedbase + 10000*(runs+1) + 1) # First "rep" of next "run"

#### Generate design matrix X for dyad covariates
Xd <- genX(n)
mX <- cbind(1,apply(Xd,3,c))
PX <- mX %*% solve(t(mX) %*% mX) %*% t(mX)     # The column space of X
# C(PAB %*% X) - the column space of X projected onto AB.
PABX <- (PAB %*% mX) %*% ginv(t(mX) %*% PAB %*% mX) %*% (t(mX) %*% PAB)
# Attempting to do a "weighted" projection
PABXw <- (PAB %*% mX) %*% diag(1 / diag(t(mX) %*% mX)) %*% (t(mX) %*% PAB)

#### Generate unobserved covariate
# ADDITIVE ROW AND COLUMN EFFECTS
true.a <- rnorm(n=n)
true.b <- rnorm(n=n)
Z.additive <- outer(true.a, true.b, "+")

Zc.low <-  eta.prime(0.0, 1, c(Z.additive), PAB, PABX)
Zc.high <- eta.prime(0.9, 1, c(Z.additive), PAB, PABX)
Z.low <-  matrix(Zc.low, nrow=n, ncol=n)
Z.high <- matrix(Zc.high, nrow=n, ncol=n)

diag(Z.additive) <- NA
diag(Z.low) <- NA
diag(Z.high) <- NA

#### Simulate binary data with Probit
eta <- intercept
for (k in 1:length(beta)) {
  eta <- eta + beta[k]*Xd[,,k]
}
eta <- eta + matrix(rnorm(n = n*n), nrow=n, ncol=n)
Y_binary_Zlow <- 1 * (eta + Z.low * gamma > 0) # Data generated with orthogonal unobserved covariate
Y_binary_Zhigh <- 1 * (eta + Z.high * gamma > 0) # Data generated with correlated unobserved covariate

# Fit a "projected" model on the orthogonal unobserved covariate
results <- amenhs::ame(Y=Y_binary_Zlow, Xdyad=Xd, family="bin",
                       R=0, rvar=TRUE, cvar=TRUE, intercept=TRUE, 
                       project=TRUE,
                       plot=FALSE, gof=FALSE, print=FALSE)

for (col in c("intercept","X1.dyad","X2.dyad","X3.dyad")) {
  png(paste(col,"lowcor_tilde.png",sep="_"), width=640, height=640)
  plot(x=results$BETA[,col], y=results$BETA.TILDE[,col],
       main=paste("Posterior Samples:",col),
       xlab="beta", ylab="beta-tilde")
  dev.off()
  print(summary(lm(results$BETA.TILDE[,col] ~ results$BETA[,col])))
  png(paste(col,"lowcor_delta.png",sep="_"), width=640, height=640)
  plot(x=results$BETA[,col], y=results$DELTA[,col],
       main=paste("Posterior Samples:",col),
       xlab="beta", ylab="delta")
  dev.off()
  print(summary(lm(results$DELTA[,col] ~ results$BETA[,col])))
}

# Fit a "projected" model on the correlated unobserved covariate
results <- amenhs::ame(Y=Y_binary_Zhigh, Xdyad=Xd, family="bin",
                       R=0, rvar=TRUE, cvar=TRUE, intercept=TRUE, 
                       project=TRUE,
                       plot=FALSE, gof=FALSE, print=FALSE)

for (col in c("intercept","X1.dyad","X2.dyad","X3.dyad")) {
  png(paste(col,"highcor_tilde.png",sep="_"), width=640, height=640)
  plot(x=results$BETA[,col], y=results$BETA.TILDE[,col],
       main=paste("Posterior Samples:",col),
       xlab="beta", ylab="beta-tilde")
  dev.off()
  print(summary(lm(results$BETA.TILDE[,col] ~ results$BETA[,col])))
  png(paste(col,"highcor_delta.png",sep="_"), width=640, height=640)
  plot(x=results$BETA[,col], y=results$DELTA[,col],
       main=paste("Posterior Samples:",col),
       xlab="beta", ylab="delta")
  dev.off()
  print(summary(lm(results$DELTA[,col] ~ results$BETA[,col])))
}

#### PLOT MEAN COMPARISONS
outdata <- as.data.frame(outdata)

for (col in c("intercept","X1.dyad","X2.dyad","X3.dyad")) {
  hc <- outdata[outdata$col==col & outdata$alpha>0.5,]
  lc <- outdata[outdata$col==col & outdata$alpha<=0.5,]
  png(paste(col,"highcor_tilde_means.png",sep="_"), width=640, height=640)
  plot(x=hc$beta.mean, y=hc$betaT.mean,
       main=paste("Posterior Means, High Correlation",col),
       xlab="beta", ylab="beta-tilde")
  abline(a=0, b=1, col="red")
  dev.off()
  png(paste(col,"highcor_delta_means.png",sep="_"), width=640, height=640)
  plot(x=hc$beta.mean, y=hc$delta.mean,
       main=paste("Posterior Means, High Correlation",col),
       xlab="beta", ylab="delta")
  abline(a=0, b=1, col="red")
  dev.off()
  png(paste(col,"lowcor_tilde_means.png",sep="_"), width=640, height=640)
  plot(x=lc$beta.mean, y=lc$betaT.mean,
       main=paste("Posterior Means, Low Correlation",col),
       xlab="beta", ylab="beta-tilde")
  abline(a=0, b=1, col="red")
  dev.off()
  png(paste(col,"lowcor_delta_means.png",sep="_"), width=640, height=640)
  plot(x=lc$beta.mean, y=lc$delta.mean,
       main=paste("Posterior Means, Low Correlation",col),
       xlab="beta", ylab="delta")
  abline(a=0, b=1, col="red")
  dev.off()
}

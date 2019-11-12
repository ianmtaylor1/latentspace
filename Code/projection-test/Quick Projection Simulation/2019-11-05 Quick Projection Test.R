library(amenhs)
library(foreach)

n <- 20       # Number of nodes in the network

runs <- 5
reps <- 10

mcmc.iter <- 10000 # Number of iterations of the MCMC algorithm to do (after "burn") (default: 10000)

intercept <- -1
beta <- c(1,1,1,1) # True coefficients (not including intercept)
gamma <- 1       # Coefficient for unobserved data

# ame() resets random seed, so seed each loop with new number
xseed <- 20191105
set.seed(xseed)

outdata <- foreach(run=1:runs, .combine="rbind") %do% {
  cat("Run:",run,"\n")
  #### Generate design matrix X for dyad covariates
  # Start empty
  Xd <- array(0, dim=c(n,n,length(beta)) )
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
  for (i in 1:n) {
    for (k in 1:length(beta)) {
      Xd[i,i,k] <- NA
    }
  }
  
  dimnames(Xd)[[3]] <- c("X1","X2","X3","X4")
  
  #### Generate unobserved covariate
  # ADDITIVE ROW AND COLUMN EFFECTS
  true.a <- rnorm(n=n)
  true.b <- rnorm(n=n)
  Z.additive <- outer(true.a, true.b, "+")
  diag(Z.additive) <- NA
  
  #### Generate unobserved covariate
  # 2-dimensional multiplicative effects
  true.u <- matrix(rnorm(n=2*n), ncol=2)
  true.v <- matrix(rnorm(n=2*n), nrow=2)
  Z.multiplicative <- true.u %*% true.v
  diag(Z.multiplicative) <- NA
  
  foreach(rep=1:reps, .combine="rbind") %do% {
    cat("Rep:",rep,"-",strftime(Sys.time(),"%H:%M:%S"),"\n")
    #### Simulate binary data with Probit
    eta <- intercept
    for (k in 1:length(beta)) {
      eta <- eta + beta[k]*Xd[,,k]
    }
    eta <- eta + matrix(rnorm(n = n*n), nrow=n, ncol=n)
    Y_binary_noZ <- 1 * (eta > 0) # Data generated WITHOUT unobserved covariate
    Y_binary_Zadditive <- 1 * (eta + Z.additive * gamma > 0) # Data generated WITH unobserved covariate
    Y_binary_Zmultiplicative <- 1 * (eta + Z.multiplicative * gamma > 0) # Data generated WITH unobserved covariate
    
    # Fit a "projected" model
    results <- amenhs::ame(Y=Y_binary_noZ, Xdyad=Xd, family="bin",
                           R=0, rvar=TRUE, cvar=TRUE, intercept=TRUE, 
                           project=TRUE,
                           plot=FALSE, gof=FALSE, print=FALSE)
    # Calculate useful things from the model
    foreach(col=c("intercept","X1.dyad","X2.dyad","X3.dyad","X4.dyad"),
            .combine="rbind") %do% {
      list(run=run,
           rep=rep,
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

write.csv(outdata, "2019-11-05 Quick Projection Results.csv", row.names=FALSE)


###############################################################################
########## JUST DO ONE TO PLOT STUFF ##########################################
###############################################################################

#### Generate design matrix X for dyad covariates
# Start empty
Xd <- array(0, dim=c(n,n,length(beta)) )
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
for (i in 1:n) {
  for (k in 1:length(beta)) {
    Xd[i,i,k] <- NA
  }
}

dimnames(Xd)[[3]] <- c("X1","X2","X3","X4")

#### Generate unobserved covariate
# ADDITIVE ROW AND COLUMN EFFECTS
true.a <- rnorm(n=n)
true.b <- rnorm(n=n)
Z.additive <- outer(true.a, true.b, "+")
diag(Z.additive) <- NA

#### Generate unobserved covariate
# 2-dimensional multiplicative effects
true.u <- matrix(rnorm(n=2*n), ncol=2)
true.v <- matrix(rnorm(n=2*n), nrow=2)
Z.multiplicative <- true.u %*% true.v
diag(Z.multiplicative) <- NA

#### Simulate binary data with Probit
eta <- intercept
for (k in 1:length(beta)) {
  eta <- eta + beta[k]*Xd[,,k]
}
eta <- eta + matrix(rnorm(n = n*n), nrow=n, ncol=n)
Y_binary_noZ <- 1 * (eta > 0) # Data generated WITHOUT unobserved covariate
Y_binary_Zadditive <- 1 * (eta + Z.additive * gamma > 0) # Data generated WITH unobserved covariate
Y_binary_Zmultiplicative <- 1 * (eta + Z.multiplicative * gamma > 0) # Data generated WITH unobserved covariate

# Fit a "projected" model
results <- amenhs::ame(Y=Y_binary_noZ, Xdyad=Xd, family="bin",
                       R=0, rvar=TRUE, cvar=TRUE, intercept=TRUE, 
                       project=TRUE,
                       plot=FALSE, gof=FALSE, print=FALSE)

for (col in c("intercept","X1.dyad","X2.dyad","X3.dyad","X4.dyad")) {
  png(paste(col,"tilde.png",sep="_"), width=640, height=640)
  plot(x=results$BETA[,col], y=results$BETA.TILDE[,col],
       main=paste("Posterior Samples:",col),
       xlab="beta", ylab="beta-tilde")
  dev.off()
  print(summary(lm(results$BETA.TILDE[,col] ~ results$BETA[,col])))
  png(paste(col,"delta.png",sep="_"), width=640, height=640)
  plot(x=results$BETA[,col], y=results$DELTA[,col],
       main=paste("Posterior Samples:",col),
       xlab="beta", ylab="delta")
  dev.off()
  print(summary(lm(results$DELTA[,col] ~ results$BETA[,col])))
}

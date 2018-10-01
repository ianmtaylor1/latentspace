library(amen)

# NOT IDEAL: I'd rather have this be in a proper library but this works for now.
source("imt_ame.R")

# Generate some test x and y
intercept <- -1
beta <- c(1,1,1) # True coefficients (not including intercept)
n <- 20

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
# Set all diagonals to NA
for (i in 1:n) {
  for (k in 1:length(beta)) {
    Xd[i,i,k] <- NA
  }
}

dimnames(Xd)[[3]] <- c("X1","X2","X3")

# Now Y
eta <- intercept
for (k in 1:length(beta)) {
  eta <- eta + beta[k]*Xd[,,k]
}
Y <- eta + matrix(rnorm(n = n*n), nrow=n, ncol=n)

# Get results!
amen.results <- ame(Y=Y, Xdyad=Xd, rvar=TRUE, cvar=TRUE, dcor=FALSE, model="nrm", gof=FALSE, plot=FALSE)
imt.results <- imt_ame(Y=Y, Xdyad=Xd)

# Compare betas
for (param in c("intercept","X1.dyad","X2.dyad","X3.dyad")) {
  par(mfcol=c(2,1))
  hist(imt.results$BETA[,param], main=paste("IMT version,",param))
  hist(amen.results$BETA[,param], main=paste("amen version,",param))
}

# Compare "Variance Components"
for (vc in c("va","vb","cab","ve")) {
  par(mfcol=c(2,1))
  hist(imt.results$VC[,vc], main=paste("IMT version, vc =",vc))
  hist(amen.results$VC[,vc], main=paste("amen version, vc =",vc))
}

# Compare A and B posterior means
imt.apm <- colMeans(imt.results$A)
imt.bpm <- colMeans(imt.results$B)
par(mfcol=c(2,1))
plot(x=amen.results$APM, y=imt.apm, 
     main="a posterior mean", xlab="amen", ylab="IMT")
abline(a=0,b=1)
plot(x=amen.results$BPM, y=imt.bpm,
     main="b posterior mean", xlab="amen", ylab="IMT")
abline(a=0,b=1)
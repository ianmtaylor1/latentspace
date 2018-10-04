library(amen)

# NOT IDEAL: I'd rather have this be in a proper library but this works for now.
source("imt_ame.R")

# Generate some test x and y
intercept <- -1
beta <- c(1,1,1) # True coefficients (not including intercept)
sigma.squared <- 1
n <- 50

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

# Generate row and column random effects
sigma.ab <- matrix(c(1,-0.5,-0.5,3),ncol=2)
ab <- rmvnorm(n=n, mu=c(0,0), Sigma=sigma.ab)
a <- ab[,1]
b <- ab[,2]

# Now Y
eta <- intercept
for (k in 1:length(beta)) {
  eta <- eta + beta[k]*Xd[,,k]
}
eta <- eta + a %*% t(rep(1,n))
eta <- eta + rep(1,n) %*% t(b)
Y <- eta + matrix(rnorm(n = n*n, mean=0, sd=sqrt(sigma.squared)), nrow=n, ncol=n)

###########
# RESULTS #
###########


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

# See how the fc updaters work with the "true" values of the parameters
fc.a <- rep(0,n)
fc.b <- rep(0,n)
fc.beta <- rep(0,length(beta)+1)
fc.sigma.squared <- 0
fc.sigma.ab <- matrix(rep(0,4), ncol=2)
# Do the reps and find the average
X <- design_array(Xdyad=Xd, intercept=TRUE, n=n)
reps <- 10000
for (i in 1:reps) {
  fc.a <- fc.a + update_a_fc(Y=Y, Z=Y, X=X, beta=c(intercept,beta), sigma.squared=sigma.squared, b=b, sigma.ab=sigma.ab)
  fc.b <- fc.b + update_b_fc(Y=Y, Z=Y, X=X, beta=c(intercept,beta), sigma.squared=sigma.squared, a=a, sigma.ab=sigma.ab)
  fc.beta <- fc.beta + update_beta_fc(Y=Y, Z=Y, X=X, sigma.squared=sigma.squared, a=a, b=b, sigma.ab=sigma.ab)
  fc.sigma.squared <- fc.sigma.squared + update_sigma.squared_fc(Y=Y, Z=Y, X=X, beta=c(intercept,beta), a=a, b=b, sigma.ab=sigma.ab)
  fc.sigma.ab <- fc.sigma.ab + update_sigma.ab_fc(Y=Y, Z=Y, X=X, beta=c(intercept,beta), sigma.squared=sigma.squared, a=a, b=b)
}
fc.a <- fc.a/reps
fc.b <- fc.b/reps
fc.beta <- fc.beta/reps
fc.sigma.squared <- fc.sigma.squared/reps
fc.sigma.ab <- fc.sigma.ab/reps
# Output the comparisons
cat("BETA: posterior mean =",fc.beta,"true =",c(intercept,beta),"\n")
cat("SIGMA.SQUARED: posterior mean =",fc.sigma.squared,"true =",sigma.squared,"\n")
cat("SIGMA.AB: posterior mean =",fc.sigma.ab,"true =",sigma.ab,"\n")
par(mfcol=c(2,1))
plot(x=a, y=fc.a, 
     main="a posterior mean", xlab="True values", ylab="FC values")
abline(a=0,b=1)
plot(x=b, y=fc.b,
     main="b posterior mean", xlab="True values", ylab="FC values")
abline(a=0,b=1)
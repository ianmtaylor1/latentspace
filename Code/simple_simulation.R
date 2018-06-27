library(amen)

n <- 20  # Size of network

intercept <- -1
beta <- c(1,1,1) # True coefficients (not including intercept)
gamma <- 1       # Coefficient for unobserved data

# Generate design matrix X for dyad covariates
# Start empty
Xd <- array(0, dim=c(n,n,length(beta)) )
for (i in 1:n) {
  for (k in 1:length(beta)) {
    Xd[i,i,k] <- NA
  }
}
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

# Generate unobserved covariate
# Similar to the shared class but with three classes (b/c why not)
uo_groups <- sample(0:1, n, TRUE) * 2 - 1
if (all(uo_groups == -1) || all(uo_groups == 1)) {
  k <- sample(1:n, 1)[1]
  uo_groups[k] <- -uo_groups[k]
}
Z <- outer(uo_groups, uo_groups, "*")
for (i in 1:n) {
  Z[i,i] <- NA
}

# Simulate binary Y data with Probit
eta <- intercept
for (k in 1:length(beta)) {
  eta <- eta + beta[k]*Xd[,,k]
}
eta <- eta + Z * gamma
eta <- eta + matrix(rnorm(n = n*n), nrow=n, ncol=n)
Y <- 1 * (eta > 0)

# Fit a model to just X, with a 1-dimensional latent multiplicative effect
fit_ame <- ame(Y, Xd=Xd, model="bin", R=1)
summary(fit_ame)


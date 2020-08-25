
library(amenhs)

set.seed(12345)

n <- 20
noise <- 0.5
intercept <- -1
beta.r <- c(1,1,1,1,1)
beta.c <- 1

# Generate row and column covariates and some 
Xr <- matrix(rnorm(5*n), ncol=5)
Xc <- rnorm(n)
a.true <- (1 - noise) * c(Xr %*% beta.r) + rnorm(n, sd=sqrt(sum(beta.r^2))*noise)
b.true <- beta.c * ((1 - noise) * Xc + rnorm(n, sd=noise))

# Response variable
Y <- (intercept + outer(c(Xr %*% beta.r), Xc * beta.c, "+") + outer(a.true, b.true, "+") + matrix(rnorm(n*n), nrow=n) >= 0)

res <- ame(Y, Xrow=Xr, Xcol=Xc, family="bin", halfcauchy=TRUE, project=TRUE,
           nscan=100000, burn=5000, odens=50,
           print=TRUE, plot=FALSE, gof=FALSE)
library(foreach)
library(batch)

jobname <- ""
n <- 20
reps <- 100
run <- 1

parseCommandArgs()

cors <- c(0)
beta <- 1
gamma <- 1

run.rnd.base <- 100000000
cor.rnd.base <- 100000
rep.rnd.base <- 1
cor.rnd.offset <- 500000

filename <- paste("latentspace_proj_corr_",jobname,"_",format(Sys.time(),"%Y-%m-%d_%H-%M"),".csv",sep="")

cat("**** RUN = ", run, "\n", sep="")
# Generate the hidden covariate (row covariate)
set.seed(run * run.rnd.base)
X.unobserved <- rnorm(n)
simulation.results <- foreach(rho = cors, .combine="rbind") %do% {
  cat("** rho = ", rho, "\n", sep="")
  # Correlated predictor
  set.seed(run * run.rnd.base + floor(rho * cor.rnd.base) + cor.rnd.offset)
  X <- c(amenproj::rmvnorm(1, rho*X.unobserved, (1-rho^2)*diag(n)))
  obs.rho <- cor(X, X.unobserved)
  foreach(j = 1:reps, .combine="rbind") %do% {
    cat("Rep = ", j, "\n", sep="")
    # Observations (Y)
    set.seed(run * run.rnd.base + floor(rho * cor.rnd.base) + cor.rnd.offset + j * rep.rnd.base)
    Z <- (X*beta + X.unobserved*gamma) %*% t(rep(1,n)) + matrix(rnorm(n^2), nrow=n, ncol=n)
    Y <- (Z >= 0) * 1.0
    # Fit it without projections
    fit <- amenproj::ame(Y=Y, Xrow=X, intercept=FALSE, family="bin", 
                         project=FALSE,
                         rvar=TRUE, cvar=FALSE, R=0, dcor=FALSE, 
                         plot=FALSE, gof=FALSE, print=FALSE)
    # Fit it with projections
    fitproj <- amenproj::ame(Y=Y, Xrow=X, intercept=FALSE, family="bin", 
                             project=TRUE,
                             rvar=TRUE, cvar=FALSE, R=0, dcor=FALSE, 
                             plot=FALSE, gof=FALSE, print=FALSE)
    # Output
    data.frame(run=run, rep=j, rho=rho, obs.rho=obs.rho, 
               beta=beta, gamma=gamma, 
               pctile=mean(fit$BETA <= beta), 
               pctile.proj=mean(fitproj$BETA <= beta),
               stringsAsFactors=FALSE)
  }
}

# Write results
write.csv(simulation.results, filename, row.names=FALSE)

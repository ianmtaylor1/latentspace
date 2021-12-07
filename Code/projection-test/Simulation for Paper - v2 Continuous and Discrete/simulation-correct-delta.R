# This file contains code to run one job of the simulation study for the
# restricted network regression paper. It is mean to be called via shell scripts
# on the summit computing environment.

library(foreach)
library(digest)
# Function in amenhs will be called directly via amenhs::func(...) to avoid
# any potential mixup with identically named functions in amen.

################################################################################
################################################################################

# Global parameters

result.basedir <- "results"
N <- 27

################################################################################
################################################################################

index.to.args <- function(idx) {
  stopifnot(idx >= 0, idx < 7 * 3 * 2 * 2 * 20)
  
  # Pull out options with successive mods and integer division
  
  # Run: 20 options
  run <- idx %% 20 + 1
  idx <- idx %/% 20
  
  # Response: 2 options
  response.idx <- idx %% 2 + 1
  response <- c("continuous", "binary")[response.idx]
  idx <- idx %/% 2
  
  # Number of random effects: 2 options
  num.re <- idx %% 2 + 1
  idx <- idx %/% 2
  
  # Random effect type: 3 options
  retype.idx <- idx %% 3 + 1
  re.type <- c("none", "invgamma", "halfcauchy")[retype.idx]
  idx <- idx %/% 3
  
  # Excess variation: 7 options
  # idx should now be between 0 and 6
  excessvar <- c("none", "high-large", "low-large", "ind-large", "high-small", "low-small", "ind-small")[idx + 1]
  
  list(run=run, response=response, num.re=num.re, re.type=re.type, excessvar=excessvar)
}



################################################################################
################################################################################

# Generate design matrix, etc according to the kind of run we're doing
# (everything except the error/response)

# Takes a matrix X and a vector Y and rescales Y in terms of its orthogonal components
# Until it has the desired canonical correlation with X. Assumes no centering needs to be done
# The returned vector will have the same length/magnitude as the input vector Y
gencancor <- function(X, Y, rho) {
  PX <- X %*% MASS::ginv(t(X) %*% X) %*% t(X)
  Ypar <- c(PX %*% Y)
  Yperp <- Y - Ypar
  Ynew <- (rho * sqrt(sum(Yperp^2)) * Ypar + sqrt(1-rho^2) * sqrt(sum(Ypar^2)) * Yperp)
  return(Ynew * sqrt(sum(Y^2)) / sqrt(sum(Ynew^2)))
}


deltas.df <- foreach(jobno=seq(0, 1679), .combine="rbind") %do% {
  
  arglist <- index.to.args(jobno)
  
  excessvar <- arglist$excessvar
  re.type <- arglist$re.type
  num.re <- arglist$num.re
  response <- arglist$response
  run <- arglist$run
  
  cat("Run index", jobno, ":", excessvar, re.type, num.re, response, run, "\n")
  
  # 0. Set the seed predictably
  design.seed <- strtoi(substr(digest(list(excessvar, num.re, response, run)), 1, 6), base=16)
  set.seed(design.seed)
  
  # 1. The design matrix and true beta values
  #    Always: intercept, dyad covariate, row covariate, column covariate
  Xd <- array(rnorm(N^2), dim=c(N, N, 1))
  Xc <- matrix(rnorm(N), ncol=1, nrow=N)
  Xr <- matrix(rnorm(N), ncol=1, nrow=N)
  
  # 2. Any excess variation due to unobserved factors
  #    First: get the correlation and magnitude of the excess variation
  if (excessvar == "none") {
    ev.cor <- 0
    ev.mag <- 0
  } else if (excessvar == "high-large") {
    ev.cor <- 0.9
    ev.mag <- 1
  } else if (excessvar == "low-large") {
    ev.cor <- 0.1
    ev.mag <- 1
  } else if (excessvar == "ind-large") {
    ev.cor <- 0
    ev.mag <- 1
  } else if (excessvar == "high-small") {
    ev.cor <- 0.9
    ev.mag <- 0.5
  } else if (excessvar == "low-small") {
    ev.cor <- 0.1
    ev.mag <- 0.5
  } else if (excessvar == "ind-small") {
    ev.cor <- 0
    ev.mag <- 0.5
  }
  #    Next: always generate a, the row/sender excess variation
  if (ev.mag > 0) {
    true.a <- gencancor(cbind(rep(1, N), Xr), rnorm(N, sd=ev.mag), ev.cor)
  } else {
    true.a <- rep(0, N)
  }
  #    Lastly, if 2 random effects, generate the column/receiver excess variation
  if ((num.re == 2) && (ev.mag > 0)) {
    true.b <- gencancor(cbind(rep(1, N), Xc), rnorm(N, sd=ev.mag), ev.cor)
  } else {
    true.b <- rep(0, N)
  }
  
  # The true values of beta, including the intercept
  beta.r <- beta.c <- beta.d <- 1
  intercept <- 1
  
  # Beta as a vector in the order it will show up in the samples
  beta <- c(intercept, beta.d, beta.r, beta.c)
  
  # Build X as a matrix for projecting excess variation to create true delta
  Xmat <- matrix(0, nrow=N^2, ncol=4)
  Xmat[,1] <- 1
  Xmat[,2] <- c(Xd[,,1])
  Xmat[,3] <- rep(Xr[,1], N)
  Xmat[,4] <- c(t(matrix(Xc[,1], nrow=N, ncol=N)))
  
  # Calculate true delta, as a vector and as individual components
  delta <- beta + c(solve(t(Xmat) %*% Xmat) %*% (t(Xmat) %*% c(outer(true.a, true.b, "+"))))
  delta.intercept <- delta[1]
  delta.d <- delta[2]
  delta.r <- delta[3]
  delta.c <- delta[4]
  
  # Construct data frame one row at a time
  data.frame(
    excessvar = excessvar,
    re.type = re.type,
    num.re = num.re,
    response = response,
    run = run,
    design.seed = design.seed,
    delta_int_true = delta.intercept,
    delta_row_true = delta.r,
    delta_col_true = delta.c,
    delta_dyad_true = delta.d,
    stringsAsFactors = FALSE
  )
}

write.csv(deltas.df, file.path(result.basedir, "true_deltas.csv"))


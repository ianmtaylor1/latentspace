# This file contains code to run one job of the simulation study for the
# restricted network regression paper. It is mean to be called via shell scripts
# on the summit computing environment.

library(foreach)
library(parallel)
library(doParallel)
library(digest)
# Function in amenhs will be called directly via amenhs::func(...) to avoid
# any potential mixup with identically named functions in amen.

################################################################################
################################################################################

# Global parameters

result.basedir <- "results"
cores <- 8
reps <- 200
N <- 27

iter <- 80000
burn <- 5000
thin <- 40

################################################################################
################################################################################

# Parse the command line arguments to tell what kind of run this is

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Must pass in only the array job id")
}

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

arglist <- index.to.args(as.integer(args[1]))

excessvar <- arglist$excessvar
re.type <- arglist$re.type
num.re <- arglist$num.re
response <- arglist$response
run <- arglist$run

cat("Run index", args[1], ":", excessvar, re.type, num.re, response, run, "\n")

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

# 0. Set the seed predictably
set.seed(strtoi(substr(digest(list(excessvar, re.type, num.re, response, run)), 1, 6), base=16))

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
true.a <- gencancor(cbind(rep(1, N), Xr), rnorm(N, sd=1) * ev.mag, ev.cor)
#    Lastly, if 2 random effects, generate the column/receiver excess variation
if (num.re == 2) {
  true.b <- gencancor(cbind(rep(1, N), Xc), rnorm(N, sd=1) * ev.mag, ev.cor)
} else {
  true.b <- rep(0, N)
}

beta.r <- beta.c <- beta.d <- 1
intercept <- 1

################################################################################
################################################################################

# Do the run

# Directory where all result RDSs and CSVs will be saved
savepath <- file.path(result.basedir, excessvar, re.type, 
                      paste0("num_re_", num.re), response, paste0("run", run))
dir.create(savepath, recursive=TRUE)


cl <- makeCluster(cores)
registerDoParallel(cl)

summary <- foreach(rep=seq_len(reps), .combine="rbind", .packages=c("digest")) %dopar% {
  # Set seed predictably for error
  set.seed(strtoi(substr(digest(list(excessvar, re.type, num.re, response, run, rep)), 1, 6), base=16))
  
  # Generate error
  Z <- intercept +
    outer(c(beta.r %*% Xr), c(beta.c %*% Xc), "+") + # Sender and receiver covariates
    outer(true.a, true.b, "+") + # Sender and receiver excess variation (possibly equal to zero)
    matrix(rnorm(N^2), nrow=N, ncol=N) # Noise
  for (i in seq_along(beta.d)) {
    Z <- Z + Xd[,,i] * beta.d[i] # Each dyadic covariate
  }
  if (response == "continuous") {
    Y <- Z
  } else {
    Y <- (Z >= 0) * 1
  }
  
  # Fit model
  res <- amenhs::ame(Y, Xrow = Xr, Xcol = Xc, Xdyad = Xd,
                     family = ifelse(response=="continuous", "nrm", "bin"),
                     halfcauchy = (re.type == "halfcauchy"),
                     project = TRUE,
                     rvar = (re.type != "none"),
                     cvar = ((re.type != "none") && (num.re == 2)),
                     dcor = FALSE,
                     nscan = iter, burn = burn, odens = thin,
                     print = FALSE, plot = FALSE, gof = FALSE)
  
  # Save the fit object as an RDS file
  saveRDS(res, file.path(savepath, paste0("rep", rep, ".RDS")))
  
  # Output one row of a dataframe to be combined by rbind/foreach
  data.frame(
    excessvar = excessvar,
    re.type = re.type,
    num.re = num.re,
    response = response,
    run = run,
    rep = rep,
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)

write.csv(summary, file.path(savepath, "summary.csv"), row.names=F)


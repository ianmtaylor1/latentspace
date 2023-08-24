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

scratch.basedir <- "/scratch/alpine/imtaylor@colostate.edu/restricted-network-regression/net-reg-proj"

result.basedir <- here::here("Code", "projection-test", "Simulation for Paper - v2 Continuous and Discrete", "results")
cores <- 16
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
  stopifnot(idx >= 0, idx < 7 * 3 * 2 * 2 * 100)
  
  # Pull out options with successive mods and integer division
  
  # Run: 20 options
  run <- idx %% 100 + 1
  idx <- idx %/% 100
  
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

jobno <- args[1]
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

################################################################################
################################################################################

# Do the run

# Directory where all result RDSs and CSVs will be saved
# Directory for intermediate RDS saves - one each rep
tmpsave.dir <- file.path(scratch.basedir, excessvar, re.type, 
                      paste0("num_re_", num.re), response, paste0("run", run))
# Location of final csv to output - one row per rep in this job
finalsavepath <- file.path(result.basedir, paste0("job_", jobno, ".csv"))

dir.create(tmpsave.dir, recursive=TRUE)
dir.create(result.basedir, recursive=TRUE)

if (!file.exists(finalsavepath)) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  #registerDoSEQ()
  
  summary <- foreach(rep=seq_len(reps), .combine="rbind", .packages=c("digest")) %dopar% {
    # Set seed predictably for error
    error.seed <- strtoi(substr(digest(list(excessvar, num.re, response, run, rep)), 1, 6), base=16)
    set.seed(error.seed)
    
    # Generate error
    Z <- intercept +
      outer(c(Xr %*% beta.r), c(Xc %*% beta.c), "+") + # Sender and receiver covariates
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
    
    if (!file.exists(file.path(tmpsave.dir, paste0("rep", rep, ".RDS")))) {
      # Fit model
      res <- amenhs::ame(Y, Xrow = Xr, Xcol = Xc, Xdyad = Xd,
                         family = ifelse(response=="continuous", "nrm", "bin"),
                         halfcauchy = (re.type == "halfcauchy"),
                         project = TRUE,
                         rvar = (re.type != "none"),
                         cvar = ((re.type != "none") && (num.re == 2)),
                         dcor = FALSE,
                         nscan = iter, burn = burn, odens = thin,
                         print = FALSE, plot = FALSE, gof = FALSE, 
                         seed = 2) # Different seed may help some of the errors?
      
      # We should have 4 columns in beta and delta
      stopifnot(ncol(res$BETA) == 4, ncol(res$DELTA) == 4)
      
      # Save the fit object as an RDS file
      # Save to scratch because the total amount of data will be HUGE (like 70 GB)
      saveRDS(res, file.path(tmpsave.dir, paste0("rep", rep, ".RDS")))
    } else {
      # If this run has been done already, just load it.
      res <- readRDS(file.path(tmpsave.dir, paste0("rep", rep, ".RDS")))
    }
    
    #### Process the result into a data frame that can be output and rbinded
    
    
    # Output one row of a dataframe to be combined by rbind/foreach
    ret <- data.frame(
      excessvar = excessvar,
      re.type = re.type,
      num.re = num.re,
      response = response,
      run = run,
      rep = rep,
      design.seed = design.seed,
      error.seed = error.seed,
      beta_int_true = intercept,
      beta_row_true = beta.r,
      beta_col_true = beta.c,
      beta_dyad_true = beta.d,
      delta_int_true = delta.intercept,
      delta_row_true = delta.r,
      delta_col_true = delta.c,
      delta_dyad_true = delta.d,
      stringsAsFactors = FALSE
    )
    
    varnames <- c("int", "row", "col", "dyad")
    for (i in 1:4) {
      name <- varnames[i]
      ret[[paste0("beta_", name, "_mean")]] <- mean(res$BETA[,i])
      ret[[paste0("beta_", name, "_var")]] <- var(res$BETA[,i])
      ret[[paste0("beta_", name, "_5q")]] <- quantile(res$BETA[,i], .05)
      ret[[paste0("beta_", name, "_95q")]] <- quantile(res$BETA[,i], .95)
      ret[[paste0("delta_", name, "_mean")]] <- mean(res$DELTA[,i])
      ret[[paste0("delta_", name, "_var")]] <- var(res$DELTA[,i])
      ret[[paste0("delta_", name, "_5q")]] <- quantile(res$DELTA[,i], .05)
      ret[[paste0("delta_", name, "_95q")]] <- quantile(res$DELTA[,i], .95)
    }
    
    # Output the data frame row
    ret
  }
  
  stopCluster(cl)
  
  write.csv(summary, finalsavepath, row.names=F)
  
} else {
  cat("Result file already exists, stopping.")
}


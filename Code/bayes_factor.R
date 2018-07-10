library(amen)

# Computes the bayes factor for comparing a model with no random effects vs
# a model with a row random effect.
# Only works for normal models, without an intercept, non-symmetric
bayes_factor <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL,
                         seed=1, samples=10000) {
  # Baked in parameters for this version
  model <- "nrm"
  intercept <- FALSE
  prior <- list()
  
  set.seed(seed)
  
  # g-prior setting for nromal data 
  if(family=="nrm" & is.null(prior$g))
  { 
    prior$g<-sum(!is.na(Y))*var(c(Y),na.rm=TRUE)
  }
  # Note: There are a bunch of prior hyperparameters that I don't think
  # are necessary for the normal model...
  
  # Build the design matrix
  X<-design_array(Xrow,Xcol,Xdyad,intercept,n)

  likelihoods_no_rvar <- c()
  likelihoods_rvar <- c()
  for (s in 1:samples) {
    # Draw sigma-squared from prior
    # Draw sigma-squared-a from prior
    # Draw a from conditional prior, given sigma-squared-a
    # Draw beta from conditional prior, given sigma-squared
    # Calculate likelihood with a
    # Calculate likelihood without a
  }
  
  return mean(likelihoods_rvar)/mean(likelihoods_no_rvar)
}


##### Main section of code

n <- 20
a_var <- 1
e_var <- 1
runs <- 10
reps <- 100

rndseedbase <- 10000000

for (run in 1:runs) {
  set.seed(rndseedbase*run)
  # Generate a
  a <- rnorm(n, 0, a_var)
  for (rep in 1:reps) {
    set.seed(rndseedbase*run + rep)
    # Generate errors
    e <- rnorm(n*n, 0, e_var)
    
    # Create two y's: with and without sender effects
    y_no_rvar <- e
    y_rvar <- e + t(sapply(a,function(x){rep(x,n)})) #Converts 'a' into a matrix with each element repeated for one row
    
    bf_no_rvar <- bayes_factor(y_no_rvar)
    bf_rvar <- bayes_factor(y_rvar)
  }
}



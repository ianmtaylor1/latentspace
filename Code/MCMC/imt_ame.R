library(amen)

n.iter <- 10000
save.every <- 25

update_beta_fc <- function(Y, X, sigma.squared, a, b, sigma.ab) {
  # Hyperparameters
  Sigma_beta <- diag(c(1,1)) # Need to change to n^2 * solve(t(X) %*% X)
  return(NULL)
}

update_sigma.squared_fc <- function(Y, X, beta, a, b, sigma.ab) {
  alpha_sigma <- 0.5
  beta_sigma <- 0.5
  return(NULL)
}

update_a_fc <- function(Y, X, beta, sigma.squared, b, sigma.ab) {
  return(NULL)
}

update_b_fc <- function(Y, X, beta, sigma.squared, a, sigma.ab) {
  return(NULL)
}

update_sigma.ab_fc <- function(Y, X, beta, sigma.squared, a, b) {
  V_sigma.ab <- diag(c(1,1))
  n_sigma.ab <- 3
  return(NULL)
}

imt_ame <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, intercept=TRUE
                    n.iter=10000, save.every=25) {
  # Set up data as appropriate
  stopifnot(length(dim(Y)) == 2)
  stopifnot(dim(Y)[1] == dim(Y)[2])
  n <- dim(Y)[1]
  X <- design_array(Xrow, Xcol, Xdyad, intercept, n)
  
  # Iterate through the gibbs sampler
  for (i in 1:n.iter) {
    # Update sigma^-2
    # Update Sigma_ab^-1
    # Update beta
    # Update a,b
    # Update latent y? (TBD: probit model)
    if (i %% save.every == 0) {
      # Save the current state
    }
  }
  # Return Samples
}

# Prior hyperparameters





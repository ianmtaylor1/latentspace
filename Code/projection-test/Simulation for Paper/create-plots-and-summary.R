library(ggplot2)


data.file <- "output/projections_for_paper_all.csv"

results <- read.csv(data.file, stringsAsFactors=FALSE)

# Add flag columns for whether CIs covered
vars <- c("beta_intercept", "beta_row1", "beta_row2", "beta_dyad1", "beta_dyad2",
          "delta_intercept", "delta_row1", "delta_row2", "delta_dyad1", "delta_dyad2")
for (var in vars) {
  covered <- paste0(var, "_covered")
  truth <- paste0(var, "_true")
  ci.low <- paste0(var, "_5q")
  ci.high <- paste0(var, "_95q")
  results[,covered] <- 1 * ((results[,truth] >= results[,ci.low]) & (results[,truth] <= results[,ci.high]))
}
no_re_vars <- c("beta_intercept", "beta_row1", "beta_row2", "beta_dyad1", "beta_dyad2")
for (var in no_re_vars) {
  covered <- paste0("no_re_", var, "_covered")
  truth <- paste0(var, "_true")
  ci.low <- paste0("no_re_", var, "_5q")
  ci.high <- paste0("no_re_", var, "_95q")
  results[,covered] <- 1 * ((results[,truth] >= results[,ci.low]) & (results[,truth] <= results[,ci.high]))
}
library(ggplot2)
library(dplyr)
library(foreach)


data.file <- "output/projections_for_paper_all.csv"

results <- read.csv(data.file, stringsAsFactors=FALSE)

plotdir <- "plots"

# Columns that define the attributes of a run
runcols <- list("run", "pr", "pd", "noise", "corr", "n")

# Add flag columns for whether CIs covered
coveredcols <- c() # Keep track of the names of the newly created columns here
vars <- c("beta_intercept", "beta_row1", "beta_row2", "beta_dyad1", "beta_dyad2",
          "delta_intercept", "delta_row1", "delta_row2", "delta_dyad1", "delta_dyad2")
for (var in vars) {
  covered <- paste0(var, "_covered")
  truth <- paste0(var, "_true")
  ci.low <- paste0(var, "_5q")
  ci.high <- paste0(var, "_95q")
  results[,covered] <- 1 * ((results[,truth] >= results[,ci.low]) & (results[,truth] <= results[,ci.high]))
  coveredcols <- c(coveredcols, covered)
}
no_re_vars <- c("beta_intercept", "beta_row1", "beta_row2", "beta_dyad1", "beta_dyad2")
for (var in no_re_vars) {
  covered <- paste0("no_re_", var, "_covered")
  truth <- paste0(var, "_true")
  ci.low <- paste0("no_re_", var, "_5q")
  ci.high <- paste0("no_re_", var, "_95q")
  results[,covered] <- 1 * ((results[,truth] >= results[,ci.low]) & (results[,truth] <= results[,ci.high]))
  coveredcols <- c(coveredcols, covered)
}

# Take the mean of the *_covered columns to find the CI coverage of that run
coverage_summary <- as.data.frame(
  results %>% 
    group_by(run, pr, pd, noise, corr, n) %>%
    summarize(beta_intercept_coverage = mean(beta_intercept_covered),
              beta_row1_coverage      = mean(beta_row1_covered),
              beta_row2_coverage      = mean(beta_row2_covered),
              beta_dyad1_coverage     = mean(beta_dyad1_covered),
              beta_dyad2_coverage     = mean(beta_dyad2_covered),
              delta_intercept_coverage = mean(delta_intercept_covered),
              delta_row1_coverage      = mean(delta_row1_covered),
              delta_row2_coverage      = mean(delta_row2_covered),
              delta_dyad1_coverage     = mean(delta_dyad1_covered),
              delta_dyad2_coverage     = mean(delta_dyad2_covered),
              no_re_beta_intercept_coverage = mean(no_re_beta_intercept_covered),
              no_re_beta_row1_coverage      = mean(no_re_beta_row1_covered),
              no_re_beta_row2_coverage      = mean(no_re_beta_row2_covered),
              no_re_beta_dyad1_coverage     = mean(no_re_beta_dyad1_covered),
              no_re_beta_dyad2_coverage     = mean(no_re_beta_dyad2_covered))
)

# Convert the data to long format for ggplot
longform <- foreach(var=c("intercept", "row1", "row2", "dyad1", "dyad2"), .combine="rbind") %do% {
  beta.colname <- paste0("beta_",var,"_coverage")
  delta.colname <- paste0("delta_",var,"_coverage")
  nore.colname <- paste0("no_re_beta_",var,"_coverage")
  
  data.frame(
    run = rep(coverage_summary$run, 3),
    pr = rep(coverage_summary$pr, 3),
    pd = rep(coverage_summary$pd, 3),
    noise = rep(coverage_summary$noise, 3),
    corr = rep(coverage_summary$corr, 3),
    n = rep(coverage_summary$n, 3),
    varname = var,
    fittype = c(rep("Unprojected", nrow(coverage_summary)),
                rep("Projected", nrow(coverage_summary)),
                rep("No Random Effect", nrow(coverage_summary))),
    coverage = c(coverage_summary[[beta.colname]],
                 coverage_summary[[delta.colname]],
                 coverage_summary[[nore.colname]]),
    stringsAsFactors = FALSE
  )
}
longform$fittype <- factor(longform$fittype, levels=c("No Random Effect", "Unprojected", "Projected"))
longform$varname <- factor(longform$varname, levels=c("intercept", "row1", "row2", "dyad1", "dyad2"))

# Produce plots
for (noiselvl in c(1, 2)) {
  for (corlvl in c(0.1, 0.9)) {
    g <- ggplot(
      longform %>% filter(noise==noiselvl, corr==corlvl), 
      aes(x=varname, y=coverage, color=fittype)
    ) +
      geom_boxplot() +
      geom_hline(yintercept=0.9, linetype="dotted") +
      ggtitle(paste("Noise =",noiselvl,"Correlation =", corlvl)) +
      theme(text=element_text(size=18))
    png(file.path(plotdir, paste0("noise", noiselvl, "_cor", corlvl, ".png")),
        width=720, height=600)
    print(g)
    dev.off()
  }
}

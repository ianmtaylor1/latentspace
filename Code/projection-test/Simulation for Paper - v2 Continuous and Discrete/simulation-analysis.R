library(foreach) 
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringi)


# Read in all csv's, append into one dataframe
allres <- foreach(f=list.files("results", pattern="^job.*\\.csv"), .combine="rbind") %do% {
  read.csv(file.path("results", f), stringsAsFactors = F)
}
# Overwrite correct true deltas
truedeltas <- read.csv(file.path("results", "true_deltas.csv"), stringsAsFactors = F)
seedsdone <- c()
for (i in seq_len(nrow(truedeltas))) {
  # For every row in the true deltas df, write their values to rows in allres with
  # the same design.seed. Avoid redundant work caused by duplicate design seeds
  if (!(truedeltas[i,"design.seed"] %in% seedsdone)) {
    rowidx <- allres[,"design.seed"] == truedeltas[i,"design.seed"]
    for (truecol in c("delta_int_true", "delta_row_true", "delta_col_true", "delta_dyad_true")) {
      allres[rowidx, truecol] <- truedeltas[i, truecol]
    }
    seedsdone <- c(seedsdone, truedeltas[i,"design.seed"])
  }
}
# Split the excess variation variable into correlation and magnitude
dashloc <- stri_locate(allres$excessvar, fixed="-")[,"start"]
allres$excessvarcor <- ifelse(is.na(dashloc), "ind", substr(allres$excessvar, start=1, stop=dashloc - 1))
allres$excessvarmag <- ifelse(is.na(dashloc), "none", substr(allres$excessvar, start=dashloc+1, stop=1000))
allres$excessvarcor <- factor(allres$excessvarcor, levels=c("ind", "low", "high"))
allres$excessvarmag <- factor(allres$excessvarmag, levels=c("large", "small", "none"))
# Other factors for ordering
allres$re.type <- factor(allres$re.type, levels=c("none", "invgamma", "halfcauchy"))

######## Means and Variances ###################################################


var <- "col"
resp <- "continuous"

# Plots of posterior means and variances for restricted vs non-network models
allres %>%
  filter(response == resp, num.re == 2) %>%
  pivot_wider(id_cols=c("excessvarcor", "excessvarmag", "num.re", "response", "run", "rep", "design.seed", "error.seed"), 
              names_from=re.type, values_from=paste0("delta_", var, "_mean")) %>%
  ggplot(aes(x=none, y=invgamma)) +
  geom_point(size=0.5) +
  facet_wrap(excessvarmag ~ excessvarcor) +
  geom_abline(slope=1, intercept=0) +
  coord_fixed() +
  theme(aspect.ratio = 1)

allres %>%
  filter(response == resp, num.re == 2) %>%
  pivot_wider(id_cols=c("excessvarcor", "excessvarmag", "response", "run", "rep", "design.seed", "error.seed"), 
              names_from=re.type, values_from=paste0("delta_", var, "_var")) %>%
  ggplot(aes(x=none, y=invgamma)) +
  geom_point(size=0.5) +
  facet_wrap(excessvarmag ~ excessvarcor) +
  geom_abline(slope=1, intercept=0) +
  coord_fixed() +
  theme(aspect.ratio = 1)


####### Coverage ###############################################################


# Create columns that tell whether the CI captured the true value
for (vname in c("int", "row", "col", "dyad")) {
  for (param in c("beta", "delta")) {
    truth <- paste0(param, "_", vname, "_true")
    ci.high <- paste0(param, "_", vname, "_95q")
    ci.low <- paste0(param, "_", vname, "_5q")
    capture <- paste0(param, "_", vname, "_captured")
    
    allres[,capture] <- ((allres[,ci.high] >= allres[,truth]) & (allres[,ci.low] <= allres[,truth])) * 1
  }
  
  truth <- paste0("delta_", vname, "_true")
  ci.high <- paste0("beta_", vname, "_95q")
  ci.low <- paste0("beta_", vname, "_5q")
  capture <- paste0("betadelta_", vname, "_captured")
  
  allres[,capture] <- ((allres[,ci.high] >= allres[,truth]) & (allres[,ci.low] <= allres[,truth])) * 1
}


# Plots of coverage
coverage_summary <- as.data.frame(
  allres %>% 
    group_by(excessvarcor, excessvarmag, re.type, num.re, response, run, design.seed) %>%
    summarize(beta_int_coverage  = mean(beta_int_captured),
              beta_row_coverage        = mean(beta_row_captured),
              beta_col_coverage     = mean(beta_col_captured),
              beta_dyad_coverage       = mean(beta_dyad_captured),
              delta_int_coverage = mean(delta_int_captured),
              delta_row_coverage       = mean(delta_row_captured),
              delta_col_coverage    = mean(delta_col_captured),
              delta_dyad_coverage      = mean(delta_dyad_captured),
              betadelta_int_coverage = mean(betadelta_int_captured),
              betadelta_row_coverage = mean(betadelta_row_captured),
              betadelta_col_coverage = mean(betadelta_col_captured),
              betadelta_dyad_coverage = mean(betadelta_dyad_captured))
)



coverage_summary %>% 
  filter(response == "binary", num.re == 2) %>%
  ggplot(aes(x=re.type, y=delta_row_coverage)) +
  geom_point(aes(color=re.type)) +
  facet_wrap(excessvarmag ~ excessvarcor) +
  geom_hline(yintercept=qbinom(0.95, 200, 0.9)/200) +
  geom_hline(yintercept=qbinom(0.05, 200, 0.9)/200)
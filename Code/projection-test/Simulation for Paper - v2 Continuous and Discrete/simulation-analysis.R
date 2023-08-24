library(foreach)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringi)
library(hexbin)

#resultdir <- here::here("Code", "projection-test", "Simulation for Paper - v2 Continuous and Discrete", "results")
resultdir <- file.path("Z:", "tmp", "rnr_results_tmp")

figsavedir <- here::here("Code", "projection-test", "Simulation for Paper - v2 Continuous and Discrete")

# Read in all csv's, append into one dataframe
allres <- foreach(f=list.files(resultdir, pattern="^job.*\\.csv"), .combine="rbind") %do% {
  cat(f, "\n")
  read.csv(file.path(resultdir, f), stringsAsFactors = F)
}

# Overwrite correct true deltas
truedeltas <- read.csv(file.path(resultdir, "true_deltas.csv"), stringsAsFactors = F) |>
  group_by(design.seed) |>
  reframe(delta_int_true=unique(delta_int_true), delta_row_true=unique(delta_row_true), delta_col_true=unique(delta_col_true), delta_dyad_true=unique(delta_dyad_true))
allres <- allres |> rows_update(truedeltas, by="design.seed", unmatched = "ignore")


# Split the excess variation variable into correlation and magnitude
dashloc <- stri_locate(allres$excessvar, fixed="-")[,"start"]
allres$excessvarcor <- ifelse(is.na(dashloc), "ind", substr(allres$excessvar, start=1, stop=dashloc - 1))
allres$excessvarmag <- ifelse(is.na(dashloc), "none", substr(allres$excessvar, start=dashloc+1, stop=1000))
allres$excessvarcor <- factor(allres$excessvarcor, levels=c("ind", "low", "high"))
allres$excessvarmag <- factor(allres$excessvarmag, levels=c("none", "small", "large"))
# Other factors for ordering
allres$re.type <- factor(allres$re.type, levels=c("none", "invgamma", "halfcauchy"))

######## Means and Variances ###################################################

var <- "col"
num.re <- 2
width <- 7
height <- 4

mylabeller <- function(x) {
  if (colnames(x) == "excessvarcor") {
    names <- rep("", length(x$excessvarcor))
    names <- replace(names, x$excessvarcor == "ind", "No Correlation")
    names <- replace(names, x$excessvarcor == "low", "Low Correlation")
    names <- replace(names, x$excessvarcor == "high", "High Correlation")
    return(list(excessvarcor = names))
  } else if (colnames(x) == "excessvarmag") {
    names <- rep("", length(x$excessvarmag))
    names <- replace(names, x$excessvarmag == "small", "Small Magnitude")
    names <- replace(names, x$excessvarmag == "large", "Large Mangitude")
    return(list(excessvarmag = names))
  }
}

# Create data for geom_text to create plot labels
labeltext <- data.frame(
  excessvarcor = factor(rep(c("ind", "low", "high"), 2), levels=c("ind", "low", "high")),
  excessvarmag = factor(c(rep("small", 3), rep("large", 3)), levels=c("small", "large")),
  label = paste0("G", seq(2, 7))
)

for (resp in c("continuous", "binary")) {
  
  # Plots of posterior means and variances for restricted vs non-network models
  allres %>%
    filter(response == resp, num.re == num.re, excessvarmag != "none") %>%
    pivot_wider(id_cols=c("excessvarcor", "excessvarmag", "num.re", "response", "run", "rep", "design.seed", "error.seed"), 
                names_from=re.type, values_from=paste0("delta_", var, "_mean")) %>%
    ggplot(aes(x=none, y=invgamma)) +
    #geom_point(size=0.5) +
    #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    geom_hex() +
    scale_fill_continuous(type = "viridis") +
    facet_grid(excessvarmag ~ excessvarcor, labeller = mylabeller) +
    geom_abline(slope=1, intercept=0) +
    theme(aspect.ratio = 0.75, axis.text.x = element_text(angle=30, hjust=1)) +
    theme_bw(base_family="serif") + 
    theme(legend.position = "none") +
    ggtitle("Receiver Covariate Posterior Means") +
    xlab("Non-network Model Posterior Mean") +
    ylab("Restricted Network Model Posterior Mean") + 
    geom_text(data=labeltext, mapping=aes(x=-Inf, y=Inf, label=label), hjust=-1, vjust=2 )
  ggsave(file.path(figsavedir, paste0(resp, "-mean-comparison.png")), width=width, height=height, units="in")
  
  allres %>%
    filter(response == resp, num.re == num.re, excessvar != "none") %>%
    pivot_wider(id_cols=c("excessvarcor", "excessvarmag", "response", "run", "rep", "design.seed", "error.seed"), 
                names_from=re.type, values_from=paste0("delta_", var, "_var")) %>%
    ggplot(aes(x=none, y=invgamma)) +
    #geom_point(size=0.5) +
    #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    geom_hex() +
    scale_fill_continuous(type = "viridis") +
    facet_grid(excessvarmag ~ excessvarcor, labeller = mylabeller) +
    geom_abline(slope=1, intercept=0) +
    theme(aspect.ratio = 0.75, axis.text.x = element_text(angle=30, hjust=1)) +
    theme_bw(base_family="serif") + 
    theme(legend.position = "none") +
    ggtitle("Receiver Covariate Posterior Variances") +
    xlab("Non-network Model Posterior Variance") +
    ylab("Restricted Network Model Posterior Variance") + 
    geom_text(data=labeltext, mapping=aes(x=-Inf, y=Inf, label=label), hjust=-1, vjust=2 )
  ggsave(file.path(figsavedir, paste0(resp, "-variance-comparison.png")), width=width, height=height, units="in")
}

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
coverage_summary$prior <- rep("", nrow(coverage_summary))
coverage_summary$prior <- replace(coverage_summary$prior, coverage_summary$re.type == "none", "No Random Effects")
coverage_summary$prior <- replace(coverage_summary$prior, coverage_summary$re.type == "invgamma", "Inverse-Gamma")
coverage_summary$prior <- replace(coverage_summary$prior, coverage_summary$re.type == "halfcauchy", "Half-Cauchy")
coverage_summary$prior <- factor(coverage_summary$prior, levels=c("No Random Effects", "Inverse-Gamma", "Half-Cauchy"))


width <- 7
height <- 6

coverage_summary %>% 
  filter(response == "binary", num.re == 2, excessvarmag != "none") %>%
  ggplot(aes(x=prior, y=delta_col_coverage)) +
  #geom_jitter(aes(color=prior), width=0.1, height=0, size=0.75) +
  geom_violin(aes(color=prior, fill=prior)) +
  facet_grid(excessvarmag ~ excessvarcor, labeller = mylabeller) +
  geom_hline(yintercept=0.9, alpha=0.4) +
  geom_hline(yintercept=qbinom(0.95, 200, 0.9)/200, linetype="dashed", alpha=0.4) +
  geom_hline(yintercept=qbinom(0.05, 200, 0.9)/200, linetype="dashed", alpha=0.4) +
  theme(aspect.ratio = 0.6, axis.text.x = element_text(angle=30, hjust=1)) +
  theme_bw(base_family="serif") + 
  theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  xlab("Model Random Effects") +
  ylab("90% Credible Interval Coverage") +
  ggtitle("Credible Interval Coverage in Restricted Binary Network Regression") +
  labs(color="Random Effect Prior", fill="Random Effect Prior") + 
  geom_text(data=labeltext, mapping=aes(x=Inf, y=-Inf, label=label), hjust=2, vjust=-1 )
ggsave(file.path(figsavedir, "binary-coverage.png"), width=width, height=height, units="in")

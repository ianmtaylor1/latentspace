library(foreach)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringi)
library(hexbin)

resultdir <- here::here("Code", "projection-test", "Simulation for Paper - v2 Continuous and Discrete", "results")
#resultdir <- file.path("Z:", "tmp", "rnr_results_tmp")

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

################################################################################
## Transform results data frame into something easier to use

deltamodels <- allres |> 
  select(excessvarmag, excessvarcor, re.type, num.re, response,run, rep, design.seed, error.seed, 
         (starts_with("delta") & !ends_with("true")) | ends_with("true")) |>
  rename_with(function(x) gsub("delta", "param", x), starts_with("delta") & !ends_with("true")) |>
  mutate(modelparam = "delta")

betamodels <- allres |>
  filter(re.type != "none") |>
  select(excessvarmag, excessvarcor, re.type, num.re, response,run, rep, design.seed, error.seed, 
         (starts_with("beta") & !ends_with("true")) | ends_with("true")) |>
  rename_with(function(x) gsub("beta", "param", x), starts_with("beta") & !ends_with("true")) |>
  mutate(modelparam = "beta")

longerres <- rbind(deltamodels, betamodels)
longerres <- longerres |>
  mutate(modelnum = case_when(re.type == "none" ~ "M1",
                              (modelparam == "beta") & (re.type == "invgamma") ~ "M2a",
                              (modelparam == "beta") & (re.type == "halfcauchy") ~ "M2b",
                              (modelparam == "delta") & (re.type == "invgamma") ~ "M3a",
                              (modelparam == "delta") & (re.type == "halfcauchy") ~ "M3b")
  )

longerres <- longerres |> 
  pivot_longer(param_int_mean:delta_dyad_true,
               names_to=c("param","pindex","summary"),
               names_pattern="(.*)_(.*)_(.*)",
               values_to="value") |>
  pivot_wider(names_from=c("param","summary"),
              values_from="value")

gc()

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
    xlab("No Random Effects Posterior Mean (M1)") +
    ylab("Restricted Network Model Posterior Mean (M3a)") + 
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
    xlab("No Random Effects Posterior Variance (M1)") +
    ylab("Restricted Network Model Posterior Variance (M3a)") + 
    geom_text(data=labeltext, mapping=aes(x=-Inf, y=Inf, label=label), hjust=-1, vjust=2 )
  ggsave(file.path(figsavedir, paste0(resp, "-variance-comparison.png")), width=width, height=height, units="in")
}

##### Coverage 2.0 #############################################################

# Create columns that tell whether the CI captured the true value
longerres <- longerres |>
  mutate(delta_captured = ((param_5q <= delta_true) & (param_95q >= delta_true)) * 1.0,
         beta_captured = ((param_5q <= beta_true) & (param_95q >= beta_true)) * 1.0)

coverage_summary_longer <- longerres |>
  group_by(excessvarmag, excessvarcor, re.type, num.re, response, run, design.seed,
           modelparam, modelnum, pindex) |>
  summarize(beta_coverage = mean(beta_captured),
            delta_coverage = mean(delta_captured))

modelcolors <- RColorBrewer::brewer.pal(name="Paired", n=5)
names(modelcolors) <- c("M2a", "M2b", "M3a", "M3b", "M1")

width <- 7
height <- 4

coverage.plot.addins.v2 <- function(gg) {
  gg + 
    geom_violin(aes(color=modelnum, fill=modelnum), bw=0.02) +
    #geom_jitter(aes(color=modelnum)) +
    facet_grid(excessvarmag ~ excessvarcor, labeller = mylabeller) +
    geom_hline(yintercept=0.9, alpha=0.4) +
    geom_hline(yintercept=qbinom(0.95, 200, 0.9)/200, linetype="dashed", alpha=0.4) +
    geom_hline(yintercept=qbinom(0.05, 200, 0.9)/200, linetype="dashed", alpha=0.4) +
    theme(aspect.ratio = 0.6, axis.text.x = element_text(angle=30, hjust=1)) + 
    coord_cartesian(ylim=c(0,1)) +
    theme_bw(base_family="serif") + 
    theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5, vjust=1)) +
    xlab("Model") +
    ylab("90% Credible Interval Coverage") +
    labs(color="Model", fill="Model") + 
    geom_text(data=labeltext, mapping=aes(x=Inf, y=-Inf, label=label), hjust=2, vjust=-1 ) +
    scale_color_manual(values=modelcolors) +
    scale_fill_manual(values=modelcolors)
}

coverage_summary_longer |>
  filter(response == "binary", num.re == 2, excessvarmag != "none", pindex == "col") |>
  ggplot(aes(x=modelnum, y=delta_coverage)) |>
  coverage.plot.addins.v2()
ggsave(file.path(figsavedir, "binary-coverage.png"), width=width, height=height, units="in")

coverage_summary_longer |>
  filter(response == "continuous", num.re == 2, excessvarmag != "none", pindex == "col") |>
  ggplot(aes(x=modelnum, y=delta_coverage)) |>
  coverage.plot.addins.v2()
ggsave(file.path(figsavedir, "continuous-coverage.png"), width=width, height=height, units="in")
  
################################################################################

longerres <- longerres |>
  mutate(delta_mean_error = param_mean - delta_true,
         beta_mean_error = param_mean - beta_true,
         delta_posterior_mse = delta_mean_error ^ 2 + param_var,
         beta_posterior_mse = beta_mean_error ^ 2 + param_var)

bias_summary_longer <- longerres |>
  group_by(excessvarmag, excessvarcor, re.type, num.re, response, run, design.seed,
           modelparam, modelnum, pindex) |>
  summarize(beta_bias = mean(beta_mean_error),
            delta_bias = mean(delta_mean_error),
            beta_mse = mean(beta_mean_error^2),
            delta_mse = mean(delta_mean_error^2),
            beta_mean_posterior_mse = mean(beta_posterior_mse),
            delta_mean_posterior_mse = mean(delta_posterior_mse))




bias_summary_longer %>% 
  filter(response == "binary", num.re == 2, excessvarmag != "none", pindex == "col") |>
  ggplot(aes(x=modelnum, y=abs(delta_bias))) +
  geom_violin(aes(color=modelnum, fill=modelnum)) +
  facet_grid(excessvarmag ~ excessvarcor, labeller = mylabeller)

bias_summary_longer %>% 
  filter(response == "continuous", num.re == 2, excessvarmag != "none", pindex == "col") |>
  ggplot(aes(x=modelnum, y=abs(delta_bias))) +
  geom_violin(aes(color=modelnum, fill=modelnum)) +
  facet_grid(excessvarmag ~ excessvarcor, labeller = mylabeller)
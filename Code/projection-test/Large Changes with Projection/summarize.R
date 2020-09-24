# Collect and summarize the output from the summit runs of this project
library(foreach)
library(plyr)

outputfolder <- "output"

# Combine CSV's, filling missing values/columns with NA's
allout <- foreach(f=list.files(outputfolder), .combine="rbind.fill") %do% {
  read.csv(file.path(outputfolder, f))
}

# First, see how much the posterior mean estimates of the random effects changes

# Look at the movement in row effects, limiting to cases where there was only
# 1 column covariate
out.pc1 <- allout[allout$pc==1,]
boxplot(relmove_a_proj ~ pr + corr, 
        data=aggregate(relmove_a_proj ~ pr + corr + run, out.pc1, mean),
        main="Row Random Effect Movement: pc=1")

# Look at the movement in column effects, limiting to cases where there were 5
# row covariates
out.pr5 <- allout[allout$pr == 5,]
boxplot(relmove_b_proj ~ pc + corr, 
        data=aggregate(relmove_b_proj ~ pc + corr + run, out.pr5, mean),
        main="Column Random Effect Movement: pr=5")

# Second, how much do the estimates of covariate effects change with/without
# projections?

# Look at movement in first row covariate, limiting to cases where there was only
# 1 column covariate
out.pc1$beta_row1_ciwidth <- out.pc1$beta_row1_95q - out.pc1$beta_row1_5q
out.pc1$delta_row1_ciwidth <- out.pc1$delta_row1_95q - out.pc1$delta_row1_5q
out.pc1$row1_proj_cishrink <- out.pc1$delta_row1_ciwidth / out.pc1$beta_row1_ciwidth
boxplot(row1_proj_cishrink ~ pr + corr,
        data=aggregate(row1_proj_cishrink ~ pr + corr + run, out.pc1, mean),
        ylim=c(0,1),
        main="Row CI Shrinkage: pc=1")

# Look at movement in first column covariate, limiting to cases where there were
# 5 row covariates
out.pr5$beta_col1_ciwidth <- out.pr5$beta_col1_95q - out.pr5$beta_col1_5q
out.pr5$delta_col1_ciwidth <- out.pr5$delta_col1_95q - out.pr5$delta_col1_5q
out.pr5$col1_proj_cishrink <- out.pr5$delta_col1_ciwidth / out.pr5$beta_col1_ciwidth
boxplot(col1_proj_cishrink ~ pc + corr,
        data=aggregate(col1_proj_cishrink ~ pc + corr + run, out.pr5, mean),
        ylim=c(0,1),
        main="Column CI Shrinkage: pr=5")

# Third, Why don't the CI's shrink much for 1 covariate, high correlation?
# Hypothesis: they're already pretty small b/c of the tendency for the covariate
# to be able to absorb the unobserved effect
boxplot(beta_row1_ciwidth ~ pr + corr,
        data=aggregate(beta_row1_ciwidth ~ pr + corr + run, out.pc1, mean),
        main="Row unprojected CI width: pc=1")
boxplot(beta_col1_ciwidth ~ pc + corr,
        data=aggregate(beta_col1_ciwidth ~ pc + corr + run, out.pr5, mean),
        main="Column unprojected CI width: pr=5")

# Fourth: how much do estimates of the variance parameters change?

# Look at sigma_a, with pc fixed = 1.
out.pc1$sigma_a_ciwidth <- out.pc1$sigma_a_95q <- out.pc1$sigma_a_5q
out.pc1$sigma_a_proj_ciwidth <- out.pc1$sigma_a_proj_95q <- out.pc1$sigma_a_proj_5q
out.pc1$sigma_a_cishrink <- out.pc1$sigma_a_proj_ciwidth / out.pc1$sigma_a_ciwidth
out.pc1$sigma_a_meanshrink <- out.pc1$sigma_a_proj_mean / out.pc1$sigma_a_mean
boxplot(sigma_a_meanshrink ~ pr + corr,
        data=aggregate(sigma_a_meanshrink ~ pr + corr + run, out.pc1, mean),
        main="Shrinkage of Mean estimates of sigma_a^2: pc=1",
        ylim=c(0,1))
boxplot(sigma_a_cishrink ~ pr + corr,
        data=aggregate(sigma_a_cishrink ~ pr + corr + run, out.pc1, mean),
        main="Shrinkage of CI of sigma_a^2: pc=1",
        ylim=c(0,1))
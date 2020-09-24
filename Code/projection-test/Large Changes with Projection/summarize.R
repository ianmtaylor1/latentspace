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


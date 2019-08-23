library(foreach)
library(ggplot2)

#indir <- "proj_corr_2019-08-22"
#outfile <- "proj_corr_2019-08-22.csv"
#
#combined <- foreach(f=list.files(indir), .combine="rbind") %do% {
#  dat <- read.csv(paste(indir,f, sep="/"), stringsAsFactors=FALSE)
#  exclude <- (dat[,"rho"] == 0) & (abs(dat[,"obs.rho"]) > 0.99)
#  dat[!exclude,]
#}
#
#write.csv(combined, outfile, row.names=FALSE)

datafile <- "proj_corr_2019-08-22.csv"

cl <- 0.9
lb <- (1-cl)/2
ub <- 1-lb

alldata <- read.csv(datafile, stringsAsFactors=FALSE)

runs <- unique(alldata[,"run"])
rhos <- unique(alldata[,"rho"])

results <- foreach(run=runs, .combine="rbind") %do% {
  foreach(rho=rhos, .combine="rbind") %do% {
    data <- alldata[(alldata[,"run"] == run) & (alldata[,"rho"] == rho),]
    data.frame(rho=rep(mean(data[,"obs.rho"]), 2),
               method=c("No Projection", "Projection"),
               coverage=c(mean((data[,"pctile"] <= ub) & (data[,"pctile"] >= lb)),
                          mean((data[,"pctile.proj"] <= ub) & (data[,"pctile.proj"] >= lb)))
               )
  }
}

ggplot(aes(x=rho, y=coverage), data=results) + geom_point(aes(color=method)) + geom_abline(slope=0, intercept=cl)

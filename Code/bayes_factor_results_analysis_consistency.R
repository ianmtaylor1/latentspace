# Read the CSV output from the Bayes Factor simulation and produce histograms

setwd("C:\\Users\\ianmt\\Documents\\Git\\latentspace\\Results")

results <- read.csv("Bayes_Factor_Consistency.csv")
runs <- unique(results[,"run"])

pdf(file="Bayes_Factor_Consistency.pdf",width=8.5,height=5.5)
repsperplot <- 10

for (run in runs) {
  reps <- unique(results[results[,"run"] == run , "rep"])
  numreps <- max(reps)
  numplots <- ((numreps - 1) %/% repsperplot) + 1
  ymin <- min(log(results[results[,"run"]==run,"bf_rvar"]),
              log(results[results[,"run"]==run,"bf_no_rvar"]))
  ymax <- min(max(log(results[results[,"run"]==run,"bf_rvar"]),
                  log(results[results[,"run"]==run,"bf_no_rvar"])),
              100)
  for (plot in 1:numplots) {
    plotidx <- (results[,"run"] == run) &
               (results[,"rep"] > repsperplot * (plot - 1)) &
               (results[,"rep"] <= repsperplot * plot)
    boxplot(log(bf_rvar)~rep, data=results[plotidx,], 
            main=paste("Run",run,": Rvar"), xlab="Rep", ylab="log(BF)")
    boxplot(log(bf_no_rvar)~rep, data=results[plotidx,], 
            main=paste("Run",run,": No Rvar"), xlab="Rep", ylab="log(BF)")
  }
}

dev.off()

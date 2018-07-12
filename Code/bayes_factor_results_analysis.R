# Read the CSV output from the Bayes Factor simulation and produce histograms

setwd("C:\\Users\\ianmt\\Documents\\Git\\latentspace\\Results")

results <- read.csv("Bayes_Factor.csv")
runs <- range(results[,"run"])

pdf(file="Bayes_Factor.pdf",width=8.5,height=5.5)

for (run in seq(runs[1],runs[2])) {
  xmin <- min(log(results[results[,"run"]==run,"bf_rvar"]),
              log(results[results[,"run"]==run,"bf_no_rvar"]))
  xmax <- min(max(log(results[results[,"run"]==run,"bf_rvar"]),
                  log(results[results[,"run"]==run,"bf_no_rvar"])),
              100)
  hist(log(results[results[,"run"]==run,"bf_rvar"]), freq=FALSE,
       main=paste("Histogram of Bayes Factors, Run =",run),
       xlab="log(bf)", xlim=c(xmin,xmax),
       col = rgb(1,0,0,0.25))
  hist(log(results[results[,"run"]==run,"bf_no_rvar"]), freq=FALSE,
       main=paste("Histogram of Bayes Factors, Run =",run),
       xlab="log(bf)", xlim=c(xmin,xmax),
       col = rgb(0,0,1,0.25), add=TRUE)
}

dev.off()

coverage <- function(df, basename, truevalue) {
  n <- dim(df)[1]
  covered <- sum((df[,paste(basename,"low",sep="_")] < truevalue) 
                 * (df[,paste(basename,"high",sep="_")] > truevalue))
  return(covered / n)
}

N <- rep(0,3*4*2*3)
var <- rep(0,3*4*2*3)
HiddenVar <- rep(0,3*4*2*3)
rand_effects <- rep(0,3*4*2*3)
ci_coverage <- rep(0,3*4*2*3)
analysis <- data.frame(N,var,HiddenVar,rand_effects)
i <- 1

truevalues <- c(-1,1,1,1)
names(truevalues) <- c("int","X1","X2","X3")

for (n in c(15,20,40)) {
  data <- read.csv(paste("simulation_results_n",n,".csv",sep=""))
  for (v in c("int","X1","X2","X3")) {
    for (h in c("noZ","yesZ")) {
      for (r in c("noRE","addRE","mulRE")) {
        basename <- paste(h,r,v,sep="_")
        cov <- coverage(data, basename, truevalues[v])
        analysis[i,"N"] <- n
        analysis[i,"var"] <- v
        analysis[i,"HiddenVar"] <- h
        analysis[i,"rand_effects"] <- r
        analysis[i,"ci_coverage"] <- cov
        i <- i + 1
      }
    }
  }
}

write.csv(analysis,"analysis.csv")


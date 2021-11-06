missing.start <- c()
missing.end <- c()
onrun <- F
for (i in 0:1679) {
  if (!file.exists(file.path("results", paste0("job_", i, ".csv")))) {
    if (!onrun) {
      missing.start <- c(missing.start, i)
      onrun <- T
    }
  } else {
    if (onrun) {
      missing.end <- c(missing.end, i-1)
      onrun <- F
    }
  }
}
if (onrun) {
  missing.end <- c(missing.end, 1679)
  onrun <- F
}

different <- missing.start != missing.end
text <- rep("", length(missing.start))
text[different] <- paste(missing.start[different], missing.end[different], sep="-")
text[!different] <- missing.start[!different]
cat(text, sep=",")

library(foreach) 

allres <- foreach(f=list.files("results"), .combine="rbind") %do% {
  read.csv(file.path("results", f), stringsAsFactors = F)
}

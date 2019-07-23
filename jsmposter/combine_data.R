library(foreach)

indir <- "data/jsm20"

outfile <- "JSM Results n=20.csv"

alldata <- foreach (f=list.files(indir,pattern=".*\\.csv"), .combine="rbind") %do% {
  read.csv(file.path(indir, f), stringsAsFactors=FALSE)
}

write.csv(alldata, outfile, row.names=FALSE)

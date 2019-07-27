library(foreach)

indir <- "data/jsm80"

outfile <- "data/JSM Results n=80.csv"

alldata <- foreach (f=list.files(indir,pattern=".*\\.csv"), .combine="rbind") %do% {
  read.csv(file.path(indir, f), stringsAsFactors=FALSE)
}

write.csv(alldata, outfile, row.names=FALSE)

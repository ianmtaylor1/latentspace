
library(doParallel)

registerDoParallel(4, cores=4)

data(IR90s)

gdp <- IR90s$nodevars[,2]
topgdp <- which(gdp >= sort(gdp,decreasing=TRUE)[30])
Y <- log(IR90s$dyadvars[topgdp,topgdp,2] + 1)

Xn <- IR90s$nodevars[topgdp,]
Xn[,1:2] <- log(Xn[,1:2])

Xd <- IR90s$dyadvars[topgdp,topgdp,c(1,3,4,5)]
Xd[,,3] <- log(Xd[,,3])


fits <- foreach (i = seq(1000,20000,1000)) %dopar% {
  cat(i,"\n")
  library(amen)
  ame(Y, Xd=Xd, Xc=Xn, Xr=Xn, R=2, plot=FALSE, print=FALSE, gof=FALSE, nscan=i, odens=i/400)
}

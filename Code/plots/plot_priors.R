library(RColorBrewer)

gfxfolder <- "graphics" # Change this if you want

dir.create(gfxfolder)
mypng <- function(filename, width=640, height=640, ...) {
  filename <- paste(gfxfolder,"/",filename,sep="")
  png(filename, width, height, ...)
}


# Get a good looking color palette
mycolors <- brewer.pal(8, "Set1")

# Density of the square of a Cauchy random variable
dcauchy2 <- function(x) {
    return(
      (x >= 0) * (1 / sqrt(x)) / (pi * (1 + x))
    )
}

# Density of Inverse Gamma
dinvgam <- function(x) {
  shape <- 0.5
  rate <- 0.5
  return(
    (x >= 0) * dgamma(1/x, shape=shape, rate=rate) / x^2
  )
}

# Plot Cauchy squared versus inverse Gamma
x <- seq(.01, 5, by=.001)
y.cauchy2 <- dcauchy2(x)
y.ig <- dinvgam(x)
mypng("Cauchy-squared.png")
plot(x, y.cauchy2, type="l", col=mycolors[1])
lines(x, y.ig, type="l", col=mycolors[2])
legend(x=median(x),y=max(y.cauchy2, y.ig),
       legend=c("Cauchy Squared", "Inverse Gamma"),
       col=mycolors[1:2], lty="solid", lwd=2)
dev.off()

# Marginal density of the horseshoe prior
dhorseshoe <- function(x) {
  # Numerically integrate over lambda for the given x
  lbound <- .01
  rbound <- 20
  step <- .01
  total <- 0
  for (lambda in seq(lbound,rbound,by=step)) {
    total <- total + dnorm(x=x, sd=sqrt(lambda)) * dcauchy2(lambda)
  }
  total[x==0] <- NA
  return(total * step)
}

# Marginal density of the normal-inverse gamma prior
digprior <- function(x) {
  # Numerically integrate over lambda for the given x
  lbound <- .01
  rbound <- 20
  step <- .01
  total <- 0
  for (lambda in seq(lbound,rbound,by=step)) {
    total <- total + dnorm(x=x, sd=sqrt(lambda)) * dinvgam(lambda)
  }
  return(total * step)
}

x <- seq(-5,5,by=.001)
y.hs <- dhorseshoe(x)
y.igprior <- digprior(x)
y.norm <- dnorm(x, sd=1)

mypng("Horseshoe.png")
plot(x,y.hs,type="l", col=mycolors[3],
     ylim=c(0,0.5))
points(x,y.igprior,type="l",col=mycolors[4])
points(x,y.norm,type="l",col=mycolors[5])
legend(x=1, y=0.5,
       legend=c("Horseshoe Prior", "Inverse Gamma Prior", "Standard Normal"),
       col=mycolors[3:5], lty="solid", lwd=2)
dev.off()


mypng("Quantiles-boxplot.png")
# Randomly generate a bunch of data from each prior, plot boxplots
n <- 10000
hssamples <- rnorm(n=n, sd=abs(rcauchy(n=n)))
igsamples <- rnorm(n=n, sd=sqrt(1/rgamma(n=n, shape=0.5, rate=0.5)))

ylim <- quantile(c(igsamples, hssamples), c(.05,.95))
par(mfrow=c(1,2))
boxplot(igsamples, ylim=ylim, main="Inverse-Gamma prior")
boxplot(hssamples, ylim=ylim, main="Horseshoe prior")
dev.off()

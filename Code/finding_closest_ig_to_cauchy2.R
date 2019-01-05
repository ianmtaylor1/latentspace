# Density of the Cauchy-squared distribution
dcauchy2 <- function(x) {
  (x > 0) / (pi * sqrt(x) * (1 + x))
}

# Cumulative distribution function of Cauchy-squared distribution
pcauchy2 <- function(x) {
  (x > 0) * 2 * atan(sqrt(x)) / pi
}

# Density of inverse gamma distribution (uses shape/rate parameterization)
dinvgamma <- function(x, alpha, beta) {
  filter <- (x > 0) & (alpha > 0) & (beta > 0)
  filter * ((beta ^ alpha) / gamma(alpha)) * (1 / x)^(alpha + 1) * exp(-beta / x)
}

# Cumulative distribution function of Inverse gamma distribution
pinvgamma <- function(x, alpha, beta) {
  pgamma(1/x, shape=alpha, rate=beta, lower.tail=FALSE)
}

# Kullback-Liebler divergence (version 1, relative to Cauchy-squared)
KL1 <- function(alpha, beta, subdivisions=1000) {
  tryCatch({
    integrate(function(x) {dcauchy2(x) * log(dcauchy2(x) / dinvgamma(x, alpha, beta))}, 0, Inf,
              subdivisions=subdivisions)$value
  }, 
  error = function(e) { Inf }
  )
}

# Kullback-Liebler divergence (version 2, relative to inverse gamma)
KL2 <- function(alpha, beta, subdivisions=1000) {
  tryCatch({
    integrate(function(x) {dinvgamma(x, alpha, beta) * log(dinvgamma(x, alpha, beta) / dcauchy2(x))}, 0, Inf,
              subdivisions=subdivisions)$value
  }, 
  error = function(e) { Inf }
  )
}

# L2 difference between densities
L2diff <- function(alpha, beta, subdivisions=1000) {
  tryCatch({
    integrate(function(x) {(dinvgamma(x, alpha, beta) - dcauchy2(x))^2}, 0, Inf,
              subdivisions=subdivisions)$value
  },
  error = function(e) { Inf }
  )
}

# L1 difference between densities
L1diff <- function(alpha, beta, subdivisions=1000) {
  tryCatch({
    integrate(function(x) {abs(dinvgamma(x, alpha, beta) - dcauchy2(x))}, 0, Inf,
              subdivisions=subdivisions)$value
  },
  error = function(e) { Inf }
  )
}

# Use optim to find optimal alpha and beta based on these four metrics
start = c(0.5,0.5)

# These three are always infinite. Proof?
#optim(start, function(x) {KL1(x[1], x[2])})
#optim(start, function(x) {KL2(x[1], x[2])})
#optim(start, function(x) {L2diff(x[1], x[2])})

l1opt <- optim(start, function(x) {L1diff(x[1], x[2])})$par

###############################################################################



# Finds difference between CDFs in "total variation": sup|f(x) - g(x)|
totalvar <- function(alpha, beta) {
  start <- 1
  res <- optimize(function(x) {abs(pcauchy2(x) - pinvgamma(x, alpha, beta))},
                  interval=c(0,10000), maximum=TRUE)
  res$objective
}

# Use optim to find optimal alpha and beta based on total variation
start = c(0.5, 0.5)

tvopt <- optim(start, function(x) {totalvar(x[1], x[2])})$par


###############################################################################
library(ggplot2)

# Plot all of these densities against each other.
x <- seq(.001, 10, .001)
ycauchy2 <- dcauchy2(x)
yig <- dinvgamma(x, 0.5, 0.5)
yigl1 <- dinvgamma(x, l1opt[1], l1opt[2])
yigtv <- dinvgamma(x, tvopt[1], tvopt[2])

cycauchy2 <- pcauchy2(x)
cyig <- pinvgamma(x, 0.5, 0.5)
cyigl1 <- pinvgamma(x, l1opt[1], l1opt[2])
cyigtv <- pinvgamma(x, tvopt[1], tvopt[2])

l1distname <- paste("L1-optimal IG(",round(l1opt[1],digits=3),",",round(l1opt[2],digits=3),")", sep="")
tvdistname <- paste("Total Variation-optimal IG(",round(tvopt[1],digits=3),",",round(tvopt[2],digits=3),")", sep="")

densities = data.frame(x=rep(x,4), 
                       pdf=c(ycauchy2, yig, yigl1, yigtv),
                       cdf=c(cycauchy2, cyig, cyigl1, cyigtv),
                       distribution=c(rep("Cauchy Squared",length(x)),
                                      rep("Standard IG(0.5,0.5)", length(x)),
                                      rep(l1distname, length(x)),
                                      rep(tvdistname, length(x))),
                       isnew=c(rep(FALSE,2*length(x)), rep(TRUE,2*length(x)))
                       )

ggplot(densities) +
  scale_colour_brewer(palette = "Set1") +
  geom_line(aes(x=x, y=pdf, color=distribution)) +
  coord_cartesian(xlim=c(0,5), ylim=c(0, 1.5)) +
  labs(title="Comparison of Densities", x="", y="")

ggplot(densities) +
  scale_colour_brewer(palette = "Set1") +
  geom_line(aes(x=x, y=cdf, color=distribution)) +
  coord_cartesian(xlim=c(0,10), ylim=c(0, 1)) +
  labs(title="Comparison of CDFs", x="", y="")
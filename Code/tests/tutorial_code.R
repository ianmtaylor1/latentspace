library(amen)

####################################
#### SECTION 1 #####################
####################################

data(IR90s)

gdp <- IR90s$nodevars[,2]
topgdp <- which(gdp >= sort(gdp,decreasing=TRUE)[30])
Y <- log(IR90s$dyadvars[topgdp,topgdp,2] + 1)

Y[1:5,1:5]

####################################
#### SUBSECTION 1.1 ################
####################################

### ---- ANOVA for trade data

Rowcountry <- matrix(rownames(Y),nrow(Y),ncol(Y))
Colcountry <- t(Rowcountry)

anova(lm(c(Y) ~ c(Rowcountry) + c(Colcountry)))

#### ---- comparison of countries in terms of row and column means

rmean <- rowMeans(Y, na.rm=TRUE)
cmean <- colMeans(Y, na.rm=TRUE)
muhat <- mean(Y, na.rm=TRUE)
ahat <- rmean - muhat
bhat <- cmean - muhat

# additive "exporter" effects
head(sort(ahat,decreasing=TRUE))

# additive "importer" effects
head(sort(bhat,decreasing=TRUE))

#### ---- covariance and correlation between row and column effects

cov(cbind(ahat,bhat))
cor(ahat,bhat)

#### ---- an estimate of dyadic covariance and correlation

## (how correlated is (i,j) with (j,i)?)
R <- Y - (muhat + outer(ahat,bhat,"+"))
cov(cbind(c(R),c(t(R))), use="complete")
cor(c(R), c(t(R)), use="complete")

#### ---- fitting the social relations model (SRM) using the default values of the ame() command
fit_SRM <- ame(Y)

mean(fit_SRM$BETA) # model-based overall mean estimate
apply(fit_SRM$VC[,1:3],2,mean) # model-based row/column mean covariance
mean(fit_SRM$VC[,4]) # model-based dyadic correlation

# Posterior means of row and column effects
fit_SRM$APM
fit_SRM$BPM

####################################
#### SUBSECTION 1.2 ################
####################################

#### ---- nodal covariates
dimnames(IR90s$nodevars)[[2]]

Xn <- IR90s$nodevars[topgdp,]
Xn[,1:2] <- log(Xn[,1:2])

#### ---- dyadic covariates
dimnames(IR90s$dyadvars)[[3]]

Xd <- IR90s$dyadvars[topgdp,topgdp,c(1,3,4,5)]
Xd[,,3] <- log(Xd[,,3])

#### ---- fitting the SRRM (social relations regression model)
fit_srrm <- ame(Y, Xd=Xd, Xr=Xn, Xc=Xn)
summary(fit_srrm)

#### ---- comparing to an ordinary linear regression approach
fit_rm <- ame(Y, Xd=Xd, Xc=Xn, Xr=Xn, rvar=FALSE, cvar=FALSE, dcor=FALSE)
summary(fit_rm)

####################################
#### SUBSECTION 1.3 ################
####################################

#### ---- fit previous srrm without the "interaction" effets (shared igo, and polity interaction)
fit_srrm0 <- ame(Y, Xd=Xd[,,1:2],Xc=Xn,Xr=Xn)
summary(fit_srrm0)

#### ---- fit an AME model using a rank 2 latent multiplicative factor
fit_ame2 <- ame(Y, Xd=Xd, Xc=Xn, Xr=Xn, R=2)
summary(fit_ame2)


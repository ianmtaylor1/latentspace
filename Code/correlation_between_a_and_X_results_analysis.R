setwd("C:\\Users\\ianmt\\Documents\\Git\\latentspace\\Results")

results <- read.csv("Correlation_Between_RowEffects_and_Covariates.csv")

pdf(file="Correlation_between_a_and_X_rowonly.pdf",width=8.5,height=5.5)

for (dattype in c("nrm","bin")) {
  residx <- (results[,"model"] == dattype)
  boxplot(X1cor~Run, data=results[residx,], 
          main=paste("Cor(a,X1); Data:",dattype), xlab="Run", ylab="Correlation")
  boxplot(X2cor~Run, data=results[residx,], 
          main=paste("Cor(a,X2); Data:",dattype), xlab="Run", ylab="Correlation")
  boxplot(X3cor~Run, data=results[residx,], 
          main=paste("Cor(a,X3); Data:",dattype), xlab="Run", ylab="Correlation")
  boxplot(Xtotcor~Run, data=results[residx,], 
          main=paste("Cor(a,Xtot); Data:",dattype), xlab="Run", ylab="Correlation")
}

dev.off()

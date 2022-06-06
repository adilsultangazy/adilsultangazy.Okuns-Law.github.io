#Clean your workspace
rm(list=ls())
dev.off()
#Press Ctrl+L


#Function to Fit the multiple linear regression model 
#with Covariance Matrix consistent under heteroskedasticity


Fit1<-function(Y,X,n,k){
  k1<-k+1
  XtX<-t(X)%*%X
  XtX_inv<-solve(XtX)
  beta_hat<-XtX_inv%*%t(X)%*%Y
  Y_hat<-X%*%beta_hat
  eps_hat<-Y-Y_hat
  eps2<-eps_hat*eps_hat
  A<-matrix(0,k1,k1)
  for (i in 1:n){
    A<-A+eps2[i]*X[i,]%*%t(X[i,])}
  factor<-n/(n-k1)
  Omega<-XtX_inv%*%A%*%XtX_inv
  #Omega<-factor*Omega
  #Standard Errors, t-statistics,P-values
  StdError<-rep(0,k1)
  tstat<-rep(0,k1)
  Pvalue<-rep(0,k1)
  for (j in 1:k1){
    StdError[j]<-sqrt(Omega[j,j])
    tstat[j]<-beta_hat[j]/StdError[j]
    Pvalue[j]<-2-2*pnorm(abs(tstat[j]),0,1)
  }
  char_vec<-c("beta_hat","Omega","eps_hat","StdErrors","tstat","Pvalue")
  out<-list("beta_hat"=beta_hat,"Omega"=Omega,"eps_hat"=eps_hat,
            "StdError"=StdError,"tstat"=tstat,"Pvalue"=Pvalue)
  return(out)
}

getwd()
setwd("C:/Users/Пользователь/Desktop/NU/2021 SPRING/ECON 301/ECON 301 RP")

data1 <- read.delim("ECON 301 RP TXT1.txt", header=T)

data1

Intercept<-rep(1,28)
lrgdp <- log(data1$rgdp)
lunemp <- log(data1$unem)
lunemp2<- log(data1$unem)^2
lunemlcurr <- log(data1$unem)*log(data1$curr)
lempl <- log(data1$empl)
linfl <- log(data1$infl)
linve <- log(data1$inve)
lcurr <- log(data1$curr)
migr <- data1$migr
loilp <- log(data1$oilp)

# Get Summary statistics: mean, min, max, std

summary(lrgdp)
sd(lrgdp, na.rm = TRUE)

summary(lunemp)
sd(lunemp, na.rm = TRUE)

summary(lunemlcurr)
sd(lunemlcurr, na.rm = TRUE)

summary(lunemp2)
sd(lunemp2, na.rm = TRUE)

summary(lempl)
sd(lempl, na.rm = TRUE)

summary(linfl)
sd(linfl, na.rm = TRUE)

summary(linve)
sd(linve, na.rm = TRUE)

summary(lcurr)
sd(lcurr, na.rm = TRUE)

summary(migr)
sd(migr, na.rm = TRUE)

summary(loilp)
sd(loilp, na.rm = TRUE)

#Plot the following graphs

plot(lunemp,lrgdp, 
     ylab = "The log of real GDP",
     xlab = "The log of unemployment rate",
     col = "black")

plot(lunemp2,lrgdp, 
     ylab = "The log of real GDP",
     xlab = "The log of unemployment rate squared",
     col = "black")

plot(lunemlcurr,lrgdp, 
     ylab = "The log of real GDP",
     xlab = "The interaction effect b/n change in unemployment and KZT/USD currency movement",
     col = "black")


#regressor matrix X1 for Model 1

X1<-cbind(Intercept, lunemp)
X1

slrgdp<-Fit1(lrgdp,X1,28,1)

slrgdp

slrgdp$Omega
slrgdp$beta_hat

printout<-cbind(slrgdp$beta_hat,slrgdp$StdError,slrgdp$tstat,slrgdp$Pvalue)

colnames(printout)<-c("beta_hat","StdError","t-stat","P-value")

options(scipen=7)
printout

boxplot(slrgdp$eps_hat)


#correlation matrix

data2 <-data.frame(lunemp, lempl, linfl, linve, lcurr, migr, loilp, lunemp2, lunemlcurr)
cor(data2)

# R squared for Model #1

SSE <- sum(slrgdp$eps_hat^2)
SSE
SST <- sum((lrgdp-mean(lrgdp, na.rm = TRUE))^2)
SST
rsquared <-1-(SSE/SST)
rsquared

#regressor matrix X2 for Model_2 and Model_2(quadratic and interaction effects)

X2<-cbind(Intercept, lunemp, linfl, linve, lcurr, migr, loilp, lunemlcurr)
X2


mlrgdp<-Fit1(lrgdp,X2,28,7)

mlrgdp

mlrgdp$Omega
mlrgdp$beta_hat

printout<-cbind(mlrgdp$beta_hat,mlrgdp$StdError,mlrgdp$tstat,mlrgdp$Pvalue)

colnames(printout)<-c("beta_hat","StdError","t-stat","P-value")

options(scipen=7)
printout

boxplot(mlrgdp$eps_hat)

#Plot Residuals


plot(mlrgdp$eps_hat,
     ylab = "Residuals",
     xlab = "Fitted Values",
     col = "black")

# R squared for Model #2

SSE2 <- sum(mlrgdp$eps_hat^2)
SSE2
SST2 <- sum((lrgdp-mean(lrgdp, na.rm = TRUE))^2)
SST2
rsquared2 <-1-(SSE2/SST2)
rsquared2


#Function for Testing a General Linear Hypothesis
WaldG<-function(Omega,HG,betahat,ThetaHip,df){
  Num<-HG%*%betahat-ThetaHip
  Vari<-solve(HG%*%Omega%*%t(HG))
  Wcalc<-t(Num)%*%Vari%*%Num
  Pvalue<-1-pchisq(Wcalc,df)
  out<-list(Wcalc,Pvalue)
  return(out)
}

#Testing Model Significance

k=7
HS<-matrix(0,k,k+1)
for(j in 1:k){ HS[j,j+1]=1}
HS

THip<-c(0,0,0,0,0,0,0)

TestSigM1<-WaldG(mlrgdp$Omega,HS,mlrgdp$beta_hat,THip,k)

TestSigM1[[1]]
TestSigM1[[2]]









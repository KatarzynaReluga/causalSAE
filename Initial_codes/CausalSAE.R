rm(list = ls())
#### IPW , IPW-EBLUP, IPW-MQ:
library("lme4")
library("sampling")
library("BinaryMQ")
library("MASS")
source("QRLM.R")
### IPW_EBLUP function ##################
ipw_mq <- function(y,x.y, x.c, w, area, ...) {
  
  }

ipw_eblup <- function(y,x.y, x.c, w, area, ...) {
  
  }



### save the errors #####################
ER_01 <- matrix(data = NA, nrow = NoSim , ncol = m) # Error Clasical method
ER_02 <- matrix(data = NA, nrow = NoSim , ncol = m) # Error EBLUP
ER_03 <- matrix(data = NA, nrow = NoSim , ncol = m) # Error MQ

### save each population true area effects #########
tau_true <- matrix(data = NA, nrow = NoSim , ncol = m)

### Population parameters ##############
sigma2u=3
sigma2e=6
sigma2B = 1 # variance of the treatment effects
sigma2V = 0.25 # variance of random effects for treatment status W
NoSim<-1000
m=50
ni=rep(5,m)
Ni=rep(100,m)
N=sum(Ni)
n=sum(ni)


for(j in 1:NoSim){ # creating the population data.frame and sampling at random NoSim times

set.seed(3*j)
  
tau_hat_01 <- vector(mode = "numeric", length = m) # ATE using only the sample (classical IPW)
tau_hat_02 <- vector(mode = "numeric", length = m) # ATE SAE-EBLUP method
tau_hat_03 <- vector(mode = "numeric", length = m) # ATE SAE-MQ method
  
out.u=0 #percentage of outliers in u
out.e=0 #percentace of outliers in e. I use usually 3%

# random effect for the outcome
u=rnorm(m,0,sqrt(sigma2u))
mean.u=9
uu <- u
k1 <- as.integer(out.u* m) 	# k1 = no. of outliers in u
u1 <- rnorm(k1, mean.u, sqrt(20))
uu[(m-k1+1):m] <- u1
u <- ifelse(out.u > 0, uu, u)
u=rep(u,times=Ni)
# error for the outcome
n1<-rbinom(N,1,out.e)
mean.e<-20
e <- (1-n1)*rnorm(N, 0, sqrt(sigma2e))+n1*rnorm(N, mean.e, sqrt(150))
gr=rep(1:m,times=Ni)
ar=unique(gr)
# random effect for treatment assignment 
V <- rnorm(m,0,sqrt(sigma2V))
V <- rep(V, times=Ni)


X=matrix(c(rlnorm(N,log(4.5)-0.5,0.5)),nrow=N,ncol=1) # extra variables that are included only in the outcome model only for precission 
Z=matrix(runif(N,min=0,max=1),nrow=N,ncol=1) # confounding factors that are included in both propensity model and outcome model.
pi<-exp(-1+0.5*Z+V)/(1+exp(-1+0.5*Z+V))
W<-rbinom(N,1,pi)
B <- rnorm(m,mean = 10, sd=sqrt(sigma2B))
B <- rep(B, times=Ni)
y=100+2*X+1*Z+ B*W +u+e # generating a scenario with different area specific effect

pop.matrix<-cbind(y,X, Z, gr,W)
pop<-as.data.frame(pop.matrix)
names(pop)<-c("y","x","z","area","w")    

# True area effects
for(oo in 1:m){
  tau_true[j,oo] <- (sum(W[pop$area==oo]*y[pop$area==oo]/pi[pop$area==oo])/sum(W[pop$area==oo]/pi[pop$area==oo]))-(sum((1-W[pop$area==oo])*y[pop$area==oo]/(1-pi[pop$area==oo]))/sum((1-W[pop$area==oo])/(1-pi[pop$area==oo])))
}

### MQ procedure ####
y_hat <- vector(mode = "numeric", length = N)
Q=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))
test="try-error"
while(test=="try-error") # It samples and fits untile it converges
{
  s <- strata(pop,"area", size=ni , method = "srswor") 
  samp <- pop[s[,2],]
  non.samp <-pop[-s[,2],]
  
  M_p <- try(glm.mq.binom(y=samp$w,x=cbind(1,samp$z),k=1.6,q=Q,maxit=100),silent=TRUE) # model to predict propensity scores
  test=class(M_p)
}

tmp.scores<-QSCORE(samp$w, M_p$fitted, Q)
scores=(cbind(samp$area,tmp.scores))
MQ_p=rep(0,m)
for (i in ar){
  MQ_p[i]=(mean(scores[,2][scores[,1]==ar[i]]))
}

M_p1 <- glm.mq.binom(y=samp$w,x=cbind(1,samp$z),k=1.6,q=MQ_p,maxit=100)

M_y <- QRLM(y=samp$y,x=cbind(1,samp$x,samp$z,samp$w),q=Q,k=1.345,maxit=100) # if we assume the treatment effect is different for each area and that the differnces are at random then we have to include the random slope for treatment status
qo<-matrix(c(gridfitinter(samp$y,M_y$fitted.values,M_y$q.values)),nrow=n,ncol=1)
qmat<-matrix(c(qo,samp$area),nrow=n,ncol=2)
Qi<-tapply(qmat[,1],qmat[,2],mean)
M_y1<-QRLM(y=samp$y,x=cbind(1,samp$x,samp$z,samp$w),q=Qi,k=1.345,maxit=100)
predict<-NULL
for (i in ar)
{ predict<-c(predict,(M_y1$coef[1,i]+M_y1$coef[2,i]*pop$x[pop$area==i]+M_y1$coef[3,i]*pop$z[pop$area==i]+M_y1$coef[4,i]*pop$w[pop$area==i]))
}

y_hat<-predict
y_hat[s[,2]] <- samp$y

p_hat<-NULL
for (i in ar)
{
  p_hat<-c(p_hat,exp(M_p1$coef[1,i]+M_p1$coef[2,i]*pop$z[pop$area==i])/(1+exp(M_p1$coef[1,i]+M_p1$coef[2,i]*pop$z[pop$area==i])))  
}


for (i in unique(pop$area)) {
  tau_hat_03[i] <- (pop$w[pop$area==i] %*% (y_hat[pop$area==i] /p_hat[pop$area==i] )) / (rep(1,Ni[i]) %*% (pop$w[pop$area==i] /p_hat[pop$area==i]))-
    ((1-pop$w[pop$area==i]) %*% (y_hat[pop$area==i] / (1-p_hat[pop$area==i]))) / (rep(1,Ni[i]) %*% ( (1-pop$w[pop$area==i]) / (1-p_hat[pop$area==i])))
}


### EBLUP procedure ####
y_hat <- vector(mode = "numeric", length = N)

M_y <- lmer(y ~ x + z + w + (1+w||area) ,data = samp) # # model to predict the y's, which leads to EBLUP
M_p <- glmer(w ~ z + (1|area) ,data = samp, family = binomial(logit)) #  model to predict propensity scores

p_hat<-NULL
p_hat <- predict(M_p, newdata=pop, type = "response", allow.new.levels=TRUE)
y_hat[s[,2]] <- samp$y
y_hat[-s[,2]] <- predict(M_y, newdata=non.samp, allow.new.levels=TRUE)

for (i in unique(pop$area)) {

# in the classical version, we cannot estimate tau_hat_01 for areas where they are all treated or all not treated (result in NaN),
    tau_hat_01[i] <- (samp$w[samp$area==i] %*% (samp$y[samp$area==i] /p_hat[s[s$area==i,2]] )) / (rep(1,ni[i]) %*% (samp$w[samp$area==i] /p_hat[s[s$area==i,2]]))-
    ((1-samp$w[samp$area==i]) %*% (samp$y[samp$area==i] / (1-p_hat[s[s$area==i,2]]))) / (rep(1,ni[i]) %*% ( (1-samp$w[samp$area==i]) / (1-p_hat[s[s$area==i,2]])))
    
    tau_hat_02[i] <- (pop$w[pop$area==i] %*% (y_hat[pop$area==i] /p_hat[pop$area==i] )) / (rep(1,Ni[i]) %*% (pop$w[pop$area==i] /p_hat[pop$area==i]))-
     ((1-pop$w[pop$area==i]) %*% (y_hat[pop$area==i] / (1-p_hat[pop$area==i]))) / (rep(1,Ni[i]) %*% ( (1-pop$w[pop$area==i]) / (1-p_hat[pop$area==i])))
}


ER_01[j, ] <- (tau_hat_01 - tau_true[j,])
ER_02[j, ] <- (tau_hat_02 - tau_true[j,])
ER_03[j, ] <- (tau_hat_03 - tau_true[j,])
print(c("01-b:",j))
}  

RB_01 <- colMeans(ER_01, na.rm = TRUE)/colMeans(tau_true)
RRMSE_01 <- sqrt(colMeans((ER_01)^2, na.rm = TRUE))/colMeans(tau_true)
RB_02 <- colMeans(ER_02)/colMeans(tau_true)
RRMSE_02 <- sqrt(colMeans((ER_02)^2))/colMeans(tau_true)
RB_03 <- colMeans(ER_03)/colMeans(tau_true)
RRMSE_03 <- sqrt(colMeans((ER_03)^2))/colMeans(tau_true)

save(tau_true, ER_01, ER_02, ER_03, RB_01,RB_02,RB_03,RRMSE_01,RRMSE_02,RRMSE_03,file = "Result_01_b.RData")
Direct <- ER_01+tau_true
EBLUP<- ER_02+tau_true
MQ <- ER_03+tau_true


for(oo in 1:50){
plot(density(rnorm(length(MQ[,oo]),10,1)), col="red", main = c("area",oo))
lines(density(MQ[,oo]),col=4,lty=4) #blue
lines(density(tau_true[,oo])) #black
lines(density(Direct[,oo],na.rm = T),col=2,lty=2) #red
lines(density(EBLUP[,oo]), col=3, lty=3)#green
}

## Analytical and block bootstrap estimation of MSE for MQ based estimator ###


rm(list = ls())
#### SETAREH:
library("lme4")
library("mbest")
library("sampling")
library("BinaryMQ")
library("MASS")
source("QRLM.R")

#Population generation
sigma2u=3
sigma2e=6
sigma2B = 1 # variance of the treatment effects
sigma2V = 0.25 # variance of random effects for treatment status W
NoSim<-1000
NoBoot <- 200
m=50
ni=rep(5,m)
Ni=rep(100,m)
N=sum(Ni)
n=sum(ni)

ER <- matrix(data = NA, nrow = NoSim , ncol = m) # Error MQ
ER_boot <- array(data=NA,dim = c(NoBoot,NoSim,m) )
stored.est_tau_boot <- array(data=NA,dim = c(NoBoot,NoSim,m) )
stored.true_tau_boot <- array(data=NA,dim = c(NoBoot,NoSim,m) )
M1 <- matrix(data=NA,nrow = NoSim, ncol = m) # FIRST PART OF VARIANCE
M2 <-  matrix(data=NA,nrow = NoSim, ncol = m) # Second PART OF VARIANCE
M3 <-  matrix(data=NA,nrow = NoSim, ncol = m) # Bias
MSE_est_MQ <- matrix(data=NA,nrow = NoSim, ncol = m)
MSE_boot <- matrix(data=NA,nrow = NoSim, ncol = m)

for(j in 1:NoSim){ # sampling at random NoSim times
  set.seed(3*j)
  
  out.u=0 #outliers in u of out.u=1
  out.e=0 #percentace of outliers in e. I use usually 3%
  #Asymmetric outliers in u
  u=rnorm(m,0,sqrt(sigma2u))
  mean.u=9
  uu <- u
  k1 <- as.integer(0.2 * m) 	# k1 = no. of outliers in u
  u1 <- rnorm(k1, mean.u, sqrt(20))
  uu[(m-k1+1):m] <- u1
  out.u <- rep(out.u, m)
  u <- ifelse(out.u > 0, uu, u)
  u=rep(u,each=100)
  #Asymmetric outliers in e
  n1<-rbinom(N,1,out.e)
  mean.e<-20
  e <- (1-n1)*rnorm(N, 0, sqrt(sigma2e))+n1*rnorm(N, mean.e, sqrt(150))
  gr=rep(1:m,each=100)
  ar=unique(gr)
  # random effect for treatment assignment 
  V <- rnorm(m,0,sqrt(sigma2V))
  V <- rep(V, each=100)
  
  
  X=matrix(c(rlnorm(N,log(4.5)-0.5,0.5)),nrow=N,ncol=1)
  Z=matrix(runif(N,min=0,max=1),nrow=N,ncol=1)
  pi<-exp(-1+0.5*Z+V)/(1+exp(-1+0.5*Z+V))
  W<-rbinom(N,1,pi)
  B <- rnorm(m,mean = 10, sd=sqrt(sigma2B))
  # tau_true <- B
  B <- rep(B, each=100)
  # y=100+2*X+1*Z+10*W+u+e
  y=100+2*X+1*Z+ B*W +u+e # generating a scenario with different area specific effect
  
  #plot(X,y)
  #plot(W,y)
  
  pop.matrix<-cbind(y,X, Z, gr,W)
  pop<-as.data.frame(pop.matrix)
  names(pop)<-c("y","x","z","area","w")    
  
  tau_true <- vector(mode = "numeric", length = m)
  for(oo in 1:m){
    tau_true[oo] <- (sum(W[pop$area==oo]*y[pop$area==oo]/pi[pop$area==oo])/sum(W[pop$area==oo]/pi[pop$area==oo]))-(sum((1-W[pop$area==oo])*y[pop$area==oo]/(1-pi[pop$area==oo]))/sum((1-W[pop$area==oo])/(1-pi[pop$area==oo])))
  }
  
  
   
  test="try-error"
  while(test=="try-error")
  {
    s <- strata(pop,"area", size=ni , method = "srswor") 
    samp <- pop[s[,2],]
    non.samp <-pop[-s[,2],]
    
    y_hat <- vector(mode = "numeric", length = N)
    tau_hat <- vector(mode = "numeric", length = m) # ATE SAE-MQ method
    Q<-c(0.25,0.30,0.4,seq(from=0.45,to=0.55,by=0.005),0.60,0.65,0.75)
    
    M_p <- try(glm.mq.binom(y=samp$w,x=cbind(1,samp$z),k=1.6,q=Q,maxit=100),silent=TRUE) # model to predict propensityscores
    test=class(M_p)
  }
  tmp.scores<-QSCORE(samp$w, M_p$fitted, Q)
  scores=(cbind(samp$area,tmp.scores))
  
  MQ_p=rep(0,m)
  for (i in ar){
    MQ_p[i]=(mean(scores[,2][scores[,1]==ar[i]]))
  }
  
  M_p1 <- glm.mq.binom(y=samp$w,x=cbind(1,samp$z),k=1.6,q=MQ_p,maxit=100)
  
  #M_p <- glmer(w ~ z + (1|area) ,data = samp, family = binomial(logit))
  #p_hat<-predict(M_p, newdata=pop, type = "response", allow.new.levels=TRUE)
  p_hat<-NULL
  u_hat.glm<-NULL
  for (i in ar){
    p_hat<-c(p_hat,exp(M_p1$coef[1,i]+M_p1$coef[2,i]*pop$z[pop$area==i])/(1+exp(M_p1$coef[1,i]+M_p1$coef[2,i]*pop$z[pop$area==i])))  
    z.U<-cbind(1,mean(pop$z[pop$area==i]))
    u_hat.glm[i]<-as.numeric(z.U%*%(M_p1$coef[,i]-M_p$coef[,which(Q==0.5)]))
  }
  # 
  
  ### MQ procedure
  Q=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))
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
  
  
  
  for (i in unique(pop$area)) {
    tau_hat[i] <- (pop$w[pop$area==i] %*% (y_hat[pop$area==i] /p_hat[pop$area==i] )) / (rep(1,Ni[i]) %*% (pop$w[pop$area==i] /p_hat[pop$area==i]))-
      ((1-pop$w[pop$area==i]) %*% (y_hat[pop$area==i] / (1-p_hat[pop$area==i]))) / (rep(1,Ni[i]) %*% ( (1-pop$w[pop$area==i]) / (1-p_hat[pop$area==i])))
  }
  
  
  ER[j, ] <- (tau_hat - tau_true)
  
  ###end################################################################################################################### 
  
  ### BlockBootstrap MSE
  y.fit <- (M_y$coef[1,which(Q==0.5)]+M_y$coef[2,which(Q==0.5)]*samp$x+M_y$coef[3,which(Q==0.5)]*samp$z+M_y$coef[4,which(Q==0.5)]*samp$w)
  gamma_hat <- vector(mode = "numeric", length = m) 
  r_hat <- samp$y - y.fit
  u_hat <- rep(0,m)
  
  for(i in 1:m){
    gamma_hat[i] <- sum((samp$w[samp$area==i]-mean(samp$w[samp$area==i]))*(r_hat[samp$area==i]-mean(r_hat[samp$area==i])))/sum(samp$w[samp$area==i]-mean(samp$w[samp$area==i])^2)
    u_hat[i]  <- mean(r_hat[samp$area==i])-gamma_hat[i]*mean(samp$w[samp$area==i])
    if(gamma_hat[i]=="NaN"){
      gamma_hat[i]=0
      u_hat[i]  <- mean(r_hat[samp$area==i])
    }
  }
  e_hat <- r_hat - samp$w*rep(gamma_hat,times=ni) - rep(u_hat,times=ni)
  
  MM <- mhglm(y ~ x + z + w + (1+w|area) , gaussian, data = samp)
  M_sigma2e <- sigma(MM)^2       # estimated sigma for error terms from the model
  M_sigma2u <- VarCorr(MM)$area[1,1]       # estimated sigma for random intercept from the model
  M_sigma2v <- VarCorr(MM)$area[2,2]
  
  #rescaling and centering
  gamma_hat <- sqrt(M_sigma2v)*gamma_hat/sd(gamma_hat)
  #y.resid.u <- psi.Huber(x=l.y.resid.u,c=BBB.c*sqrt(l.estsigma2u.m),MAD=FALSE)
  gamma_hat<- gamma_hat-mean(gamma_hat)
  
  u_hat <- sqrt(M_sigma2u)*u_hat/sd(u_hat)
  #y.resid.u <- psi.Huber(x=l.y.resid.u,c=BBB.c*sqrt(l.estsigma2u.m),MAD=FALSE)
  u_hat <- u_hat-mean(u_hat)
  
  e_hat <- sqrt(M_sigma2e)*e_hat/sd(e_hat)
  #y.resid.e <- psi.Huber(x=l.y.resid.e,c=BBB.c*sqrt(l.estsigma2e.m),MAD=FALSE)
  e_hat <- e_hat-mean(e_hat)
  
  GLMM <- mhglm( w ~1+ z +(1|area) , binomial(link = "logit"), data = pop)
  GM_sigma2u <- VarCorr(GLMM)$area[1,1]
  u_hat.glm<-sqrt(GM_sigma2u)*u_hat.glm/sd(u_hat.glm)
  u_hat.glm <- u_hat.glm-mean(u_hat.glm)
  
  for(i in 1:NoBoot)
  {
    set.seed(46*i)
    y_hat_boot <- vector(mode = "numeric", length = N)
    tau_boot <- vector(mode = "numeric", length = m)
    tau_true_boot <- vector(mode = "numeric", length = m)
    id.gamma<- sample(1:m, size=m ,replace = TRUE)
    gamma_boot <- (gamma_hat[id.gamma])
    u_hat.glm.boot<-u_hat.glm[id.gamma]
    u_boot <- rep(sample(u_hat, size=m ,replace = TRUE), times=Ni)
    e_boot <- vector(mode = "numeric")
    p_hat_boot <- vector(mode = "numeric")
    for(oo in 1:m){
      #id.Z<- sample(1:m, size=1 ,replace = TRUE)
      hh<-exp(M_p$coef[1,14]+M_p$coef[2,14]*pop$z[pop$area==oo]+u_hat.glm.boot[oo])
      p_hat_boot <- c(p_hat_boot,(hh/(1+hh)))
    }
    Wboot<-rbinom(N,1,p_hat_boot)
    
    for(oo in 1:m){
      h <- sample(1:m, 1)
      hh <- sample(e_hat[samp$area==h], size=Ni[oo] , replace=TRUE)
      e_boot <- c(e_boot,hh)
    }
    
    y_boot <- NULL
    for (oo in 1:m){
      y_boot <- c(y_boot,M_y$coef[1,which(Q==0.5)]+M_y$coef[2,which(Q==0.5)]*pop$x[pop$area==oo]+M_y$coef[3,which(Q==0.5)]*pop$z[pop$area==oo]+(M_y$coef[4,which(Q==0.5)]+gamma_boot[oo])*Wboot[pop$area==oo]+u_boot[pop$area==oo]+e_boot[pop$area==oo])
    }
    
    boot.data <- cbind.data.frame(y_boot , xboot=pop$x , zboot=pop$z , wboot=Wboot , areaboot=pop$area)
    Q.boot<-c(0.25,0.30,0.4,seq(from=0.45,to=0.55,by=0.005),0.60,0.65,0.75)
    
    
    test="try-error"
    test1="try-error"
    while(test=="try-error" |test1=="try-error" ){
      s_boot <- strata(boot.data,"areaboot", size=ni , method = "srswor") 
      samp.boot <- boot.data[s_boot[,2],]
      non.samp.boot <-boot.data[-s_boot[,2],]
      M_p_boot  <- try(glm.mq.binom(y=samp.boot$wboot,x=cbind(1,samp.boot$zboot),k=1.6,q=Q.boot,maxit=100),silent=TRUE) # model to predict propensityscores
      test=class(M_p_boot)
      
      tmp.scores<-QSCORE(samp.boot$wboot, M_p_boot$fitted, Q.boot)
      scores=(cbind(samp.boot$areaboot,tmp.scores))
      
      MQ_p_boot=rep(0,m)
      for (oo in 1:m){
        MQ_p_boot[oo]=(mean(scores[,2][scores[,1]==ar[oo]]))
      }
      
      M_p1_boot <- try(glm.mq.binom(y=samp.boot$wboot,x=cbind(1,samp.boot$zboot),k=1.6,q=MQ_p_boot,maxit=100),silent=TRUE)
      test1=class(M_p1_boot)
    }
    # 
    
    #M_p_boot <- glmer(wboot ~ zboot + (1|areaboot) ,data = samp.boot, family = binomial(logit))
    #p_hat_boot<-predict(M_p_boot, newdata=boot.data, type = "response", allow.new.levels=TRUE)
    est.p_hat_boot <- NULL
    for (oo in 1:m)
    {
      est.p_hat_boot<-c(est.p_hat_boot,exp(M_p1_boot$coef[1,oo]+M_p1_boot$coef[2,oo]*boot.data$zboot[boot.data$areaboot==oo])/(1+exp(M_p1_boot$coef[1,oo]+M_p1_boot$coef[2,oo]*boot.data$zboot[boot.data$areaboot==oo])))  
    }  
    
    Q=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))
    M_y_boot <- QRLM(y=samp.boot$y_boot,x=cbind(1,samp.boot$xboot,samp.boot$zboot,samp.boot$wboot),q=Q,k=1.345,maxit=100) # if we assume the treatment effect is different for each area and that the differnces are at random then we have to include the random slope for treatment status
    qo<-matrix(c(gridfitinter(samp.boot$y_boot,M_y_boot$fitted.values,M_y_boot$q.values)),nrow=n,ncol=1)
    qmat<-matrix(c(qo,samp.boot$areaboot),nrow=n,ncol=2)
    Qi<-tapply(qmat[,1],qmat[,2],mean)
    M_y1_boot<-QRLM(y=samp.boot$y_boot,x=cbind(1,samp.boot$xboot,samp.boot$zboot,samp.boot$wboot),q=Qi,k=1.345,maxit=100)
    predict_boot<-NULL
    for (oo in ar){ 
      predict_boot<-c(predict_boot,(M_y1_boot$coef[1,oo]+M_y1_boot$coef[2,oo]*boot.data$xboot[boot.data$areaboot==oo]+M_y1_boot$coef[3,oo]*boot.data$zboot[boot.data$areaboot==oo]+M_y1_boot$coef[4,oo]*boot.data$wboot[boot.data$areaboot==oo]))
    }
    y_hat_boot <- predict_boot
    y_hat_boot[s_boot[,2]] <- samp.boot$y_boot
    for (ii in 1:m) {
      tau_boot[ii] <- (boot.data$wboot[boot.data$areaboot==ii] %*% (y_hat_boot[boot.data$areaboot==ii] /est.p_hat_boot[boot.data$areaboot==ii] )) / (rep(1,Ni[ii]) %*% (boot.data$wboot[boot.data$areaboot==ii] /est.p_hat_boot[boot.data$areaboot==ii]))-
        ((1-boot.data$wboot[boot.data$areaboot==ii]) %*% (y_hat_boot[boot.data$areaboot==ii] / (1-est.p_hat_boot[boot.data$areaboot==ii]))) / (rep(1,Ni[ii]) %*% ( (1-boot.data$wboot[boot.data$areaboot==ii]) / (1-est.p_hat_boot[boot.data$areaboot==ii])))
    }
    
    for (ii in 1:m) {
      tau_true_boot[ii] <- (boot.data$wboot[boot.data$areaboot==ii] %*% (boot.data$y_boot[boot.data$areaboot==ii] /p_hat_boot[boot.data$areaboot==ii] )) / (rep(1,Ni[ii]) %*% (boot.data$wboot[boot.data$areaboot==ii] /p_hat_boot[boot.data$areaboot==ii]))-
        ((1-boot.data$wboot[boot.data$areaboot==ii]) %*% (boot.data$y_boot[boot.data$areaboot==ii] / (1-p_hat_boot[boot.data$areaboot==ii]))) / (rep(1,Ni[ii]) %*% ( (1-boot.data$wboot[boot.data$areaboot==ii]) / (1-p_hat_boot[boot.data$areaboot==ii])))
    }
    
    ER_boot[i,j,] <- tau_true_boot -tau_boot
    
    stored.est_tau_boot[i,j,] <-tau_boot
    stored.true_tau_boot[i,j,] <-tau_true_boot
    
  }
  
  MSE_boot[j,] <- apply(ER_boot[,j,]^2,2,mean,na.rm=TRUE)
  print(c("01_b:",j))
}  



save(ER , ER_boot, MSE_est_MQ,MSE_boot, file = "MSE_MQ_01_b_double_block_bootstrap_test_nicola.RData")

MSE.true <- sqrt(apply(ER^2, 2, mean,na.rm=TRUE))
MSE.est <- sqrt(apply((MSE_boot), 2, mean,na.rm=TRUE))
summary(((MSE.est-MSE.true)/MSE.true)*100)
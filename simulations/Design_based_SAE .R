###################################################################################################################
#### 
###################################################################################################################
# Changes compare to previous version : instead of using the income we consider the log of equivalized income.
# Age is restricted to be between 25 and 80 included.
# I added the family type to the Propensity score equation
# The common support within each area is considered
# The balance of the covariates via propensity scores is checked
## The following warning is repeated many times: boundary (singular) fit: see ?isSingular 
###################################################################################################################

library("lme4")
library("sampling")
library("BinaryMQ")
library("lattice")
library("GLMMadaptive")
library("xtable")
library("MatchIt")
library("ggplot2")

rm(list=ls())
load("data_eusilc15.RData")
source("QRLM.R")

y=log(data$HS130) # Lowest monthly income to make end meet
familysize=data$HX050 # Equivalized household size
y=y/familysize # log of Equivalized monthly income to make end meet

treatment=data$PL140-1# Treatment dummy where 0: Permanent , 1:temporary
age=data$RX010 # Age later trimmed to be above 24 years
gender=factor(data$RB090) #1: Male , 2: Female
education=factor(data$PE040) # only 6 categories (there are three missing categories?)
statciv = factor(data$PB190) # marital status
familytype=factor(data$TF)
#italian= factor(data$NCITT) #being the italian citizen or not (is omitted due to many missings)
tenure= factor(data$HX070) #having a tenure or not 
htype = factor(data$HX060) #household type which is not the same as family type
noroom= data$HH030
dtype= factor(data$HH010)
crime= factor(data$HS190)
dincome=data$HY020
hincome=data$fylav
wintens = data$RX040 
Xind <- cbind(age, gender, education, statciv, tenure,wintens) 
Xhh <- cbind(familytype, noroom,dtype,crime, dincome) 

region<-data$REGIONE
regname <- as.numeric(names(table(region)))
m <- length(regname)
for(i in 1:m){
  region[region==regname[i]] <- i 
}
## we put Bolzano and Trento together 
## we put Abruzzo and Molise

region[region==5] <- 4
region[region==15] <- 14

###################################################################################
## omit the missings from the dataframe
###################################################################################
dd <- data.frame(y=y, Xind, Xhh, treatment, region)
dd <- dd[order(dd$region),]
### table of missing values ####
colSums(is.na(dd)) ## first treatment( 5559), then italian(803) , third dtype(61)

dd <- na.omit(dd)
dd <- dd[dd$age>24,] # we loose 2810
## we omit Friuli-Venezia-Giulia and Liguria 
## these are special cases, that is region==c(7,8)
## dd <- dd[!(dd$region==7 | dd$region==8), ]
ar<-sort(unique(dd$region)) #id region
###################################################################################
## standardizing the continuios variables 
###################################################################################
dd$age <- scale(dd$age, center = T , scale = T)
dd$wintens <- scale(dd$wintens, center = T, scale = T)
dd$noroom <- scale(dd$noroom,center = T,scale = T)
dd$dincome <-scale(dd$dincome,center = T,scale = T)


####################################################################################
## First fit for the propensity scores
####################################################################################
modglme <- glmer(treatment ~ age+factor(gender)+factor(education)+factor(statciv)+
                  factor(tenure)+factor(familytype)+(1|region),
                data = dd, family = binomial(link="logit"),nAGQ = 0)
p_hat<-NULL
p_hat <- predict(modglme, newdata=dd, type = "response", allow.new.levels=TRUE)
summary(p_hat[dd$treatment==1])
summary(p_hat[dd$treatment==0])
 
####################################################################################
## The Common Support Assumption within each area
####################################################################################
ind <- NULL
   
 for(i in ar){
 maxi <- min(max(p_hat[dd$treatment==0 & dd$region==i]),max(p_hat[dd$treatment==1 & dd$region==i]))
 mini <- max(min(p_hat[dd$treatment==1 & dd$region==i]),min(p_hat[dd$treatment==0 & dd$region==i]))
 ind <- c(ind, (p_hat[dd$region==i]<=maxi)*(p_hat[dd$region==i]>=mini))
 }

p_hat <- subset(p_hat,ind==1 )
dd<- subset(dd,ind==1 )
 
par(mai=c(1.5,1.5,1.5,1.5))
pdf(file="Common_support.pdf", width=10)
plot(density(p_hat[dd$treatment==0]),
     col="red",  type = "l", lwd=2, main="Common support in the full sample after trimming",
     xlab="Propensity scores",ylab="Density", cex.lab=1.5,cex.main=1.5,lty="dashed", font.lab=2)
lines(density(p_hat[dd$treatment==1]), col="blue", type = "l", lwd=2)
legend(0.15,8.5,legend = c("treated","control"),lty=c(1,2), bty ="n", col=c("blue", "red"), cex=1.5, lwd = c(2,2))
dev.off()
#####################################################################################

#####################################################################################
## Full sample (population) size
#####################################################################################
Ni<-as.numeric(table(dd$region))
Nc <- as.numeric(table(dd$region[dd$treatment==0]))
Nt <- as.numeric(table(dd$region[dd$treatment==1]))
N<-sum(Ni)

#####################################################################################
## Assesing the balance in covariates
#####################################################################################
l.p.hat <- log(p_hat/(1-p_hat))
l.c.bar <- tapply(l.p.hat[dd$treatment==0], dd$region[dd$treatment==0], mean) 
l.t.bar <- tapply(l.p.hat[dd$treatment==1], dd$region[dd$treatment==1], mean)
s.l.c <- tapply(l.p.hat[dd$treatment==0], dd$region[dd$treatment==0], sd) 
s.l.t <- tapply(l.p.hat[dd$treatment==1], dd$region[dd$treatment==1], sd)

delta.l= (l.t.bar-l.c.bar)/sqrt((s.l.c^2+s.l.t^2)/2)
nu1 <- (s.l.c^2/Nc+s.l.t^2/Nt)^2
nu2 <- (s.l.c^2/Nc)^2/(Nc-1) + (s.l.t^2/Nt)^2/(Nt-1)
nu <- nu1/nu2
p.val <- 2*(1-pt(delta.l,df=nu))
t.crit <- qt(0.975,df=nu)
delta.l < t.crit  
## TRUE means that the balance is covariates is respected 8there no diffrence in 
## the mean of the log scores within each area
####################################################################################

##################################################################################
## Assessing the balance in the covariates by matchit package to reply to the referrees
##################################################################################
for(i in ar){
sub_dd <- dd[dd$region==i,]
s.out <- matchit(treatment ~ age+factor(gender)+factor(education)+factor(statciv)+
                   factor(tenure)+factor(familytype),"subclass", subclass = 10 ,data = sub_dd)
print(i)
print(summary(s.out))
}

s.out <- matchit(treatment ~ age+factor(gender)+factor(education)+factor(statciv)+
                   factor(tenure)+factor(familytype),"subclass", subclass = 10 ,data = dd)
summary(s.out)

s.out <- matchit(treatment ~ age+factor(gender)+factor(education)+factor(statciv)+
                   factor(tenure)+factor(familytype) ,"subclass", subclass = 10 ,data = dd)
summary(s.out)
####################################################################################
## Ture values of treatment 
####################################################################################
True<-NULL
 for (i in ar) {
   True[i] <- (dd$treatment[dd$region==i] %*% (dd$y[dd$region==i] /p_hat[dd$region==i] )) / (sum(dd$treatment[dd$region==i] /p_hat[dd$region==i]))-
   ((1-dd$treatment[dd$region==i]) %*% (dd$y[dd$region==i] / (1-p_hat[dd$region==i]))) / (sum( (1-dd$treatment[dd$region==i]) / (1-p_hat[dd$region==i])))
 }
True <- na.omit(True)
 
True_total <- (dd$treatment %*% (dd$y /p_hat )) / (sum(dd$treatment /p_hat))-
   ((1-dd$treatment) %*% (dd$y / (1-p_hat))) / (sum( (1-dd$treatment) / (1-p_hat)))

## Original estimates of the effects
tau.original <- rbind(cbind(Ni,Nc,Nt,na.omit(True)), c(N,sum(Nc),sum(Nt),True_total))
NUT2 <- c("Piemonte","Valle-d'Aosta","Lombardia","Bolzano-Trento","Veneto","Friuli-Venezia-Giulia", "Liguria" ,
         "Emilia-Romagna", "Toscana", "Umbria","Marche","Lazio",
         "Abruzzo-Molise","Campania","Puglia", "Basilicata",
         "Calabria","Sicilia","Sardegna", "Italy") 
est.orig <- data.frame( "NUTS2"=NUT2, "Original sample size"=tau.original[,1], 
         "Treated percentage"=round(tau.original[,3]/tau.original[,2],2), "Original estimate"=round(tau.original[,4],2))
org.table <- xtable(est.orig)
align(org.table) <- xalign(org.table)
digits(org.table) <- xdigits(org.table)
display(org.table) <- xdisplay(org.table)
print(org.table,include.rownames = FALSE)

## to show that the true values of the impact are far from being homogenious 
pdf(file="hist_true_log.pdf", width=10)
par(mai=c(1.5,1.5,1.5,1.5))
hist(True,freq = F , xlab=expression(tau[j]), cex.lab=1.5,font.lab=2, main =NULL)
lines((density(True)),col="red",lty="dashed", lwd=2)
dev.off()
####################################################################################


####################################################################################
## Design based simulations
####################################################################################
NoSim<-1000
EBLUP<-matrix(NA,NoSim,length(ar))
MQ<-matrix(NA,NoSim,length(ar))
Direct<-matrix(NA,NoSim,length(ar))
mse_Direct <- matrix(NA,NoSim,length(ar))

ni<-as.integer(0.10*Ni)
n<-sum(ni)

#### Starting simulation #####
for (h in 1:NoSim)
{
  
  tau_hat_01 <- vector(mode = "numeric", length = length(ar)) # ATE using only the sample (classical approach)
  tau_hat_02 <- vector(mode = "numeric", length = length(ar)) # ATE SAE-EBLUP method
  tau_hat_03 <- vector(mode = "numeric", length = length(ar)) # ATE SAE-MQ method
  
  set.seed(h*3)
  # draw a sample and apply the methodology
  test <- "try-error"
  test1 <- "try-error"
  test2 <- 0
  s <-NULL
    while(test == "try-error" | test1 == "try-error"|test2 <32)
  {
    s <- strata(dd,"region", size=ni , method = "srswor") 
    samp <- dd[s[,2],]
    non.samp <-dd[-s[,2],]
    ## MQ-IPW
    Q<-c(0.25,0.30,0.4,seq(from=0.45,to=0.55,by=0.005),0.60,0.65,0.75)
    X <- cbind(1,samp$age,model.matrix(~ factor(samp$gender))[,-1], model.matrix(~ factor(samp$education))[,-1], model.matrix(~factor(samp$statciv))[,-1], model.matrix(~factor(samp$tenure))[,-1],model.matrix(~ factor(samp$familytype))[,-1]) 
    pX <- cbind(1,dd$age,model.matrix(~ factor(dd$gender))[,-1], model.matrix(~ factor(dd$education))[,-1], model.matrix(~factor(dd$statciv))[,-1], model.matrix(~factor(dd$tenure))[,-1],model.matrix(~ factor(dd$familytype))[,-1]) 
    m_p <- try(glm.mq.binom(y=dd$treatment,x=pX,k=1.6,q=Q,maxit=100),silent=TRUE)
    test <- class(m_p)
    tmp.scores<-try(QSCORE(dd$treatment, m_p$fitted, Q),silent=TRUE)
    scores=try((cbind(dd$region,tmp.scores)),silent=TRUE)
    MQ_p=rep(0,length(ar))
    for (i in ar){
      MQ_p[which(ar==i)]=try((mean(scores[,2][scores[,1]==i])),silent=TRUE)
    }
    m_p1 <- try(glm.mq.binom(y=dd$treatment,x=pX,k=1.6,q=MQ_p,maxit=100),silent=TRUE)
    test1 <- class(m_p1)
    Xy <-cbind(1,samp$age,model.matrix(~ factor(samp$gender))[,-1], model.matrix(~ factor(samp$education))[,-1], model.matrix(~factor(samp$statciv))[,-1], model.matrix(~factor(samp$tenure))[,-1], samp$wintens, model.matrix(~factor(samp$familytype))[,-1] , samp$noroom ,model.matrix(~factor(samp$dtype))[,-1], model.matrix(~factor(samp$crime))[,-1], samp$dincome, samp$treatment) 
    test2 <- ncol(Xy)
  }
  
  pXy <- cbind(1,dd$age,model.matrix(~ factor(dd$gender))[,-1], model.matrix(~ factor(dd$education))[,-1], model.matrix(~factor(dd$statciv))[,-1], model.matrix(~factor(dd$tenure))[,-1], dd$wintens, model.matrix(~factor(dd$familytype))[,-1] , dd$noroom ,model.matrix(~factor(dd$dtype))[,-1], model.matrix(~factor(dd$crime))[,-1], dd$dincome, dd$treatment) 
  Q=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))
  m_y <- QRLM(y=samp$y,x=Xy,q=Q,k=1.345,maxit=100) # if we assume the treatment effect is different for each area and that the differnces are at random then we have to include the random slope for treatment status
  qo<-matrix(c(gridfitinter(samp$y,m_y$fitted.values,m_y$q.values)),nrow=n,ncol=1)
  qmat<-matrix(c(qo,samp$region),nrow=n,ncol=2)
  Qi<-tapply(qmat[,1],qmat[,2],mean)
  m_y1<-QRLM(y=samp$y,x=Xy,q=Qi,k=1.345,maxit=100)
  
  predict<-NULL
  for (i in 1:length(ar))
  { predict<-c(predict,t(m_y1$coef[,i])%*%t(pXy[dd$region==ar[i],]))
  }
  
  y_hat<-predict
  y_hat[s[,2]] <- samp$y

  p_hat<-NULL
   for (i in 1:length(ar))
   {
     index <- t(m_p1$coefficients[,i])%*%t(pX[dd$region==ar[i],])
     p_hat<-c(p_hat , exp(index)/(1+exp(index)))
   }

  
  for (i in 1:length(ar)) {
    tau_hat_03[i] <- (dd$treatment[dd$region==ar[i]] %*% (y_hat[dd$region==ar[i]] /p_hat[dd$region==ar[i]] )) / (sum(dd$treatment[dd$region==ar[i]] /p_hat[dd$region==ar[i]]))-
      ((1-dd$treatment[dd$region==ar[i]]) %*% (y_hat[dd$region==ar[i]] / (1-p_hat[dd$region==ar[i]]))) / (sum( (1-dd$treatment[dd$region==ar[i]]) / (1-p_hat[dd$region==ar[i]])))
    
  }
  
  MQ[h,]<-tau_hat_03
  
  
  ## Classic-IPW and EBLUP-IPW
  y_hat <- vector(mode = "numeric", length = N)
  
  modlme <- lmer(y~ age +factor(gender)+factor(education)+factor(statciv)+factor(tenure)+ wintens+ factor(familytype)+ noroom+factor(dtype)+factor(crime)+dincome+treatment+(1+treatment||region), data = samp)
  modglme <- glmer(treatment ~ age+factor(gender)+factor(education)+factor(statciv)+factor(tenure)+ factor(familytype) +(1|region), data = dd, family = binomial(logit),nAGQ = 0)
  
  p_hat<-NULL
  p_hat <- predict(modglme, newdata=dd, type = "response", allow.new.levels=TRUE)
  y_hat[s[,2]] <- samp$y
  y_hat[-s[,2]] <- predict(modlme, newdata=non.samp, allow.new.levels=TRUE)
  
  for (i in 1:length(ar)) {
    tau_hat_01[i] <- (samp$treatment[samp$region==ar[i]] %*% (samp$y[samp$region==ar[i]] /p_hat[s[s$region==ar[i],2]] )) / (sum(samp$treatment[samp$region==ar[i]] /p_hat[s[s$region==ar[i],2]]))-
      ((1-samp$treatment[samp$region==ar[i]]) %*% (samp$y[samp$region==ar[i]] / (1-p_hat[s[s$region==ar[i],2]]))) / (sum( (1-samp$treatment[samp$region==ar[i]]) / (1-p_hat[s[s$region==ar[i],2]])))
    
    tau_hat_02[i] <- (dd$treatment[dd$region==ar[i]] %*% (y_hat[dd$region==ar[i]] /p_hat[dd$region==ar[i]] )) / (sum(dd$treatment[dd$region==ar[i] ] /p_hat[dd$region==ar[i]]))-
      ((1-dd$treatment[dd$region==ar[i]]) %*% (y_hat[dd$region==ar[i]] / (1-p_hat[dd$region==ar[i]]))) / (sum( (1-dd$treatment[dd$region==ar[i]]) / (1-p_hat[dd$region==ar[i]])))
  }
  
  Direct[h,]<-tau_hat_01
  EBLUP[h,]<-tau_hat_02
  print(c(h,":", sum(is.na(tau_hat_01)) ,":", sum((table(samp$region,by=samp$treatment)==0)) ))
} #simulation ends


save(NUT2,EBLUP,MQ,Direct,True,dd, p_hat, file="Eusilc_simulation_03_log(y).RData")

cbind( "NUTS2"=NUT2[-20], True,"IPW_Direct"=apply(Direct,2,mean,na.rm=TRUE), "IPW_EBLUP"=apply(EBLUP,2,mean),"IPW_MQ"=apply(MQ,2,mean))
RB.EBLUP <- (apply(EBLUP,2,mean)-True)/abs(True)*100
RB.MQ <- (apply(MQ,2,mean)-True)/abs(True)*100
RB.Direct <- (apply(Direct,2,mean,na.rm=TRUE)-True)/abs(True)*100
summary(RB.Direct); summary(RB.EBLUP); summary(RB.MQ)

RRMSE.Direct <- (sqrt(apply(apply(Direct,1,function(x){(x-True)^2}),1,mean,na.rm=TRUE))/abs(True))*100
RRMSE.EBLUP <- (sqrt(apply(apply(EBLUP,1,function(x){(x-True)^2}),1,mean,na.rm=TRUE))/abs(True))*100
RRMSE.MQ <- (sqrt(apply(apply(MQ,1,function(x){(x-True)^2}),1,mean,na.rm=TRUE))/abs(True))*100
summary(RRMSE.Direct); summary(RRMSE.EBLUP); summary(RRMSE.MQ)

## table of estimates and their s.e.
results <- NULL
results <- data.frame( "NUTS2"=NUT2[-20], "True"=na.omit(True),"IPW_Direct"=apply(Direct,2,mean,na.rm=TRUE),"(s.e.Direct)"=sqrt(apply(apply(Direct,1,function(x){(x-True)^2}),1,mean,na.rm=TRUE)),
                         "IPW_EBLUP"=apply(EBLUP,2,mean),"(s.e.EBLUP)"=sqrt(apply(apply(EBLUP,1,function(x){(x-True)^2}),1,mean,na.rm=TRUE)),
                       "IPW_MQ"=apply(MQ,2,mean),"(s.e.MQ)"=sqrt(apply(apply(MQ,1,function(x){(x-True)^2}),1,mean,na.rm=TRUE)))
est.results <- xtable(results,digits=2)
print(est.results, include.rownames=FALSE)
##############################################################

## table of 95% confidence intervals of the estimates
results <- NULL
CI_Direct <-apply(Direct,2,function(x) quantile(x,probs = c(0.025,0.975),na.rm=TRUE))
CI_EBLUP <- apply(EBLUP,2,function(x) quantile(x,probs = c(0.025,0.975),na.rm=TRUE))
CI_MQ <-apply(MQ,2,function(x) quantile(x,probs = c(0.025,0.975),na.rm=TRUE))
results <- data.frame( "NUTS2"=NUT2[-20], "True"=na.omit(True),"LL_D"=CI_Direct[1,],
                       "UL_D"=CI_Direct[2,],"LL_E"=CI_EBLUP[1,],"UL_E"=CI_EBLUP[2,], 
                       "LL_MQ"= CI_MQ[1,], "UL_MQ" = CI_MQ[2,])
est.results <- xtable(results,digits=2)
print(est.results, include.rownames=FALSE)


pdf(file="CI_compare.pdf", width = 20, height = 10)
ggplot(data = results)+
  geom_point(mapping = aes(x=NUTS2,y=True), size=5)+
  geom_hline(yintercept=0, color="red", lty="dashed")+
  geom_errorbar(mapping = aes(x=NUTS2, ymin=LL_D, ymax=UL_D), lty="dashed")+
  geom_errorbar(mapping = aes(x=NUTS2, ymin=LL_E, ymax=UL_E), color="red")+
  geom_errorbar(mapping = aes(x=NUTS2, ymin=LL_MQ, ymax=UL_MQ), color="blue")+
  theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.y = element_text(size = rel(1.5), angle = 00))+
  ylab("Treatment effects")+
  theme(axis.title.y = element_text(size = 24))+
  theme(axis.title.x = element_text(size = 24))
dev.off()

    
###############################################################

## boxplot of the true and the mean of the estimates ####
pdf(file = "App_ATE_log.pdf", width = 20)
results < NULL
results <- data.frame("True"=True,"IPW_Direct"=apply(Direct,2,mean,na.rm=TRUE), "IPW_EBLUP"=apply(EBLUP,2,mean),"IPW_MQ"=apply(MQ,2,mean))
boxplot(results, cex.axis=1.5, font.lab=2, lwd=1.5)
dev.off()
## boxplot of RB and RRMSE : I redone this with ggplot furter down ####
pdf(file = "app_rb_rrmse_log.pdf", width=20)
par(mfrow=c(1,2))
boxplot(cbind(Direct=RB.Direct, EBLUP=RB.EBLUP,MQ=RB.MQ), main="Relative Bias", cex.axis=1.5, font.lab=2, lwd=1.5)
abline(h=0, col="red", lty="dashed",lwd=2)
boxplot(cbind(Direct=RRMSE.Direct,EBLUP=RRMSE.EBLUP,MQ=RRMSE.MQ), main="Relative Root Mean Square Error", cex.axis=1.5, font.lab=2, lwd=1.5)
abline(h=0, col="red", lty="dashed",lwd=2)
dev.off()

## Capturing the heterogeniety of the effects ####
min.x <- min(c(True,Direct,EBLUP,MQ),na.rm=T)
max.x <- max(c(True,Direct,EBLUP,MQ),na.rm=T)
max.y <- 2.5
pdf(file = "hetero-est-log.pdf", width=15)
plot(density(True),xlim=c(-1.5,1.5),ylim=c(0,max.y), lty=2, col="black", xlab = "area-specific effect",main = "Distribution of area-specific ATE - different estimators",lwd=2)
lines(density(Direct,na.rm = T), col="grey", lwd=2)
lines(density(EBLUP,na.rm=T),col="red", lwd=2)
lines(density(MQ,na.rm = T),co="blue", lwd=2)
legend(1,2.5,legend = c("True","IPW_Driect","IPW_EBLUP","IPW_MQ"), col=c("black","grey","red","blue"), lty = c(2,1,1,1),lwd=2, bty = "n")
dev.off()

# Capturing the heterogeniety of the effects
D <- apply(Direct,2,mean,na.rm=T)
E <- apply(EBLUP,2,mean)
M <- apply(MQ,2,mean)
min.x <- min(c(True,D,E,M),na.rm=T)
max.x <- max(c(True,D,E,M),na.rm=T)
max.y <- 2.5
pdf(file = "hetero-est-log2.pdf", width=15)
plot(density(True),xlim=c(-1.1,1.1),ylim=c(0,max.y), lty=2, col="black", xlab = "area-specific effect",  main = "Distribution of area-specific ATE - different estimators",lwd=2)
lines(density(D,na.rm = T), col="grey", lwd=2)
lines(density(E,na.rm=T),col="red", lwd=2)
lines(density(M,na.rm = T),co="blue", lwd=2)
legend(0.5,2.5,legend = c("True","IPW_Driect","IPW_EBLUP","IPW_MQ"), col=c("black","grey","red","blue"), lty = c(2,1,1,1),lwd=2, bty = "n")
dev.off()
## Boxplot of the effects and the summary of the results
pdf(file="VariationinOfDifferentEstimator.pdf", width = 15)
boxplot(cbind(Direct=c(Direct),EBLUP=c(EBLUP),MQ=c(MQ)),main="SATE/PATE", cex.axis=1.5, font.lab=2, lwd=1.5)
dev.off()
results <- NULL
results <- rbind(summary(apply(Direct,2,var, na.rm=T)),
                 summary(apply(EBLUP,2,var, na.rm=T)),
                 summary(apply(MQ,2,var, na.rm=T)))
colnames(results) <- c("min","Q1","Median","Mean","Q3","max") 
rownames(results) <- c("IPW-Direct","IPW-EBLUP","IPW-MQ")
est.results <- xtable(results,digits=4)
print(est.results, include.rownames=T)


## Formal kolmogrov-smirnov tes
ks.test(D,True)
ks.test(E,True)
ks.test(M,True)


## calculating the efficiency ####
RMSE.EBLUP<-sqrt(apply((EBLUP-True)^2,2,mean))
RMSE.MQ<-sqrt(apply((MQ-True)^2,2,mean))
RMSE.Direct<-sqrt(apply((Direct-True)^2,2,mean,na.rm=TRUE))

cat_ni<-NULL
cat_ni[ni<= 30]<-1
cat_ni[(ni>30 & ni<= 70)]<-2
cat_ni[ni> 70]<-3

EFF.EBLUP<-RMSE.EBLUP/RMSE.Direct*100
EFF.MQ<-RMSE.MQ/RMSE.Direct*100
results <- NULL
results <- rbind(summary(EFF.EBLUP), summary(EFF.MQ))
eff <- est.results <- xtable(results,digits=2)
print(est.results, include.rownames=FALSE)

mean(EFF.EBLUP)
mean(EFF.MQ)

### to detect the problem: ####
par(mfrow=c(3,1))
boxplot(Direct-True,use.cols = F, ylim=c(-5,5))
boxplot(EBLUP-True,use.cols = F, ylim=c(-5,5))
boxplot(MQ-True,use.cols = F, ylim=c(-5,5))

ar <- unique(dd$region)
plot(ar,True,type="b")
lines(ar, apply(Direct,2,mean,na.rm=TRUE), type="b", col="red" )
lines(ar, apply(EBLUP,2,mean), type="b", col="darkgreen")
lines(ar,apply(MQ,2,mean), type="b", col="darkblue")


## estimation of the effect for each area with the 3 method ??? ####
library(lattice) 
data <- NULL
data <- data.frame(y=c(matrix(Direct,I(1000*19),1),matrix(EBLUP,I(1000*19),1),matrix(MQ,I(1000*19),1)),reg=rep(rep(NUT2[-20],each=1000),3), meth= rep(c("Direct","EBLUP","MQ"),each=I(1000*19)))
data$reg <- factor(data$reg)#, levels = c(19,20,21,16,17,18,13,14,15,10,11,12,7,8,9,4,5,6,1,2,3))
#ind <- c(19,16,17,18,13,14,15,10,11,12,7,8,9,4,5,6,1,2,3)
bwplot(y~meth| reg, data=data, ylab=expression(bold("Causal effects")), xlab=expression(bold("IPW Estimators")), main="Causal effects evaluated for each region in Italy", layout=c(3,7),
  panel=function(...) {
  panel.abline(h=0, col="red",lty="dashed")
  panel.bwplot(...)
}) 

### Boxplots of the RB and RRMSE with ggplot####
load("~/Dropbox/IEforSAE/Application_EUSILC/Eusilc_simulation_03_log(y).RData")
ER_01 <- apply(Direct, MARGIN = 1,FUN = function(x) (x-True))
RB_01 <- rowMeans(ER_01, na.rm = T)/abs(True)*100
RRMSE_01 <- sqrt(rowMeans(ER_01^2, na.rm = T))/abs(True)*100
result <- data.frame("RB"=RB_01, "RRMSE"=RRMSE_01, "method"="Direct")
ER_02 <- apply(EBLUP, MARGIN = 1,FUN = function(x) (x-True))
RB_02 <- rowMeans(ER_02, na.rm = T)/abs(True)*100
RRMSE_02 <- sqrt(rowMeans(ER_02^2, na.rm = T))/abs(True)*100
result <- rbind(result, data.frame("RB"=RB_02, "RRMSE"=RRMSE_02, "method"="EBLUP"))
ER_03 <- apply(MQ, MARGIN = 1,FUN = function(x) (x-True))
RB_03 <- rowMeans(ER_03, na.rm = T)/abs(True)*100
RRMSE_03 <- sqrt(rowMeans(ER_03^2, na.rm = T))/abs(True)*100
result <- rbind(result, data.frame("RB"=RB_03, "RRMSE"=RRMSE_03, "method"="MQ"))

library(ggplot2)
library(cowplot)
RB_Compare <- ggplot(result, aes(x=method, y=RB), size=3)+
     geom_boxplot()+
     geom_hline(yintercept=0, color="red", lty="dashed", size=1.5)+
     ylab("Relative Bias%") + xlab("")+
     theme(text=element_text(size=24,  family="serif"))+
     theme(axis.text.x = element_text(size = rel(1.8), angle = 00))+
     theme(axis.text.y = element_text(size = rel(1.5), angle = 00))



RRMSE_Compare <- ggplot(result, aes(x=method, y=RRMSE), size=3)+
  geom_boxplot()+
  geom_hline(yintercept=0, color="red", lty="dashed", size=1.5)+
  ylab("Relative Root Mean Square Error%") + xlab("")+
  theme(text=element_text(size=24,  family="serif"))+
  theme(axis.text.x = element_text(size = rel(1.8), angle = 00))+
  theme(axis.text.y = element_text(size = rel(1.5), angle = 00))

pdf(file = "RB_RRMSE_Compare.pdf", width = 20, height = 10)
plot_grid (RB_Compare, RRMSE_Compare, ncol=2)
dev.off()

## interaction plot 
ggplot(data=dd, aes(x=treatment,y=y, group=region, color=as.factor(region)))+
  theme(axis.text.x = element_blank())+
  labs(color="Regions")+
  stat_summary(fun=mean, geom="line")

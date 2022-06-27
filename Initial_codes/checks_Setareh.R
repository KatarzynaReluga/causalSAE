set.seed(1)
#Population generation
sigma2u=3
sigma2e=6
sigma2B = 1 # variance of the treatment effects
sigma2V = 0.25 # variance of random effects for treatment status W
NoSim<-1000
NoBoot <- 200
m=50
#ni=rep(5, m)
#ni = rep(20, m)
Ni=rep(100,m)
N=sum(Ni)
n=sum(ni)


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

pop.matrix<-cbind(y,X, Z, gr,W,pi)
pop<-as.data.frame(pop.matrix)
names(pop)<-c("y","x","z","area","w", "p_score")    

tau_true1 <- vector(mode = "numeric", length = m)
tau_true2 <- vector(mode = "numeric", length = m)
for(oo in 1:m){
  tau_true1[oo] <- (sum(W[pop$area==oo]*y[pop$area==oo]/pi[pop$area==oo])/sum(W[pop$area==oo]/pi[pop$area==oo]))-(sum((1-W[pop$area==oo])*y[pop$area==oo]/(1-pi[pop$area==oo]))/sum((1-W[pop$area==oo])/(1-pi[pop$area==oo])))
  tau_true2[oo] <- mean(W[pop$area==oo]*y[pop$area==oo]/pi[pop$area==oo])-  mean((1-W[pop$area==oo])*y[pop$area==oo]/(1-pi[pop$area==oo]))
  #tau_true[oo] <- (sum(A[pop$group==oo]*y[pop$group==oo]/pi[pop$group==oo])/sum(A[pop$group==oo]/pi[pop$group==oo]))-(sum((1-A[pop$group==oo])*y[pop$group==oo]/(1-pi[pop$group==oo]))/sum((1-A[pop$group==oo])/(1-pi[pop$group==oo])))
  #tau_true[oo] <- (mean(A[pop$group==oo]*y[pop$group==oo]/pi[pop$group==oo]))-(mean((1-A[pop$group==oo])*y[pop$group==oo]/(1-pi[pop$group==oo])))
  
}
#all assumptions we make, should we rather compute means of 
tau_true1
tau_true2

pop2 <- pop
names(pop2)<-c("y","x","z","group","A", "p_score")    
A = W
true_tau <-  calculate_tau(list(pop2))[[1]]$tau
true_tau

# EBP
s <- strata(pop,"area", size=ni , method = "srswor") 
samp <- pop[s[,2],]
non.samp <-pop[-s[,2],]

y_hat <- vector(mode = "numeric", length = N)

M_y <- lmer(y ~ x + z + w + (1+w||area),data = samp) # # model to predict the y's, which leads to EBLUP
M_p <- glmer(w ~ z + (1|area) ,data = samp, family = binomial) #  model to predict propensity scores

p_hat<-NULL
p_hat <- predict(M_p, newdata=pop, type = "response", allow.new.levels=TRUE)
y_hat[s[,2]] <- samp$y
y_hat[-s[,2]] <- predict(M_y, newdata=non.samp, allow.new.levels=TRUE)

tau_hat_01 <- vector(mode = "numeric", length = m) # ATE using only the sample (classical IPW)
tau_hat_02 <- vector(mode = "numeric", length = m) # ATE SAE-EBLUP method
tau_hat_03 <- vector(mode = "numeric", length = m) # ATE SAE-MQ method

for (i in unique(pop$area)) {
  
  # in the classical version, we cannot estimate tau_hat_01 for areas where they are all treated or all not treated (result in NaN),
#  tau_hat_01[i] <- (samp$w[samp$area==i] %*% (samp$y[samp$area==i] /p_hat[s[s$area==i,2]] )) / (rep(1,ni[i]) %*% (samp$w[samp$area==i] /p_hat[s[s$area==i,2]]))-
#    ((1-samp$w[samp$area==i]) %*% (samp$y[samp$area==i] / (1-p_hat[s[s$area==i,2]]))) / (rep(1,ni[i]) %*% ( (1-samp$w[samp$area==i]) / (1-p_hat[s[s$area==i,2]])))
  
  tau_hat_01[i] <- (pop$w[pop$area==i] %*% (y_hat[pop$area==i] /p_hat[pop$area==i] )) / (rep(1,Ni[i]) %*% (pop$w[pop$area==i] /p_hat[pop$area==i]))-
    ((1-pop$w[pop$area==i]) %*% (y_hat[pop$area==i] / (1-p_hat[pop$area==i]))) / (rep(1,Ni[i]) %*% ( (1-pop$w[pop$area==i]) / (1-p_hat[pop$area==i])))

#  tau_hat_01[i] <- mean(samp$w[samp$area==i] * (samp$y[samp$area==i] /p_hat[s[s$area==i,2]] )) -
#    mean((1-samp$w[samp$area==i]) * (samp$y[samp$area==i] / (1-p_hat[s[s$area==i,2]])))
  
  tau_hat_02[i] <- mean(pop$w[pop$area==i] * (y_hat[pop$area==i] /p_hat[pop$area==i] ))-
    mean((1-pop$w[pop$area==i]) * (y_hat[pop$area==i] / (1-p_hat[pop$area==i]))) 
  
}

plot(1:m, tau_true1, type = "l")
lines(tau_hat_01, col = 2)

plot(tau_true2, type = "l")
lines(tau_hat_02, col = 3)

lines(true_tau, col = 4)
lines(ht_EBLUP$tau_hat[[1]]$tau, col = 2)
####

pop2 <- pop
names(pop2)<-c("y","x","z","group","A", "p_score")    
A = W
true_tau <-  calculate_tau(list(pop2))[[1]]$tau
true_tau

# EBP
#s2 <- strata(pop2,"group", size=ni , method = "srswor") 
s2 <- s
colnames(s2)[1] <- group
samp2 <- pop2[s2[,2],]
non.samp2 <-pop2[-s2[,2],]




ht_EBLUP <- hte(formula_y = y ~ x + z + A + (1+A||group), 
                formula_p_score = A ~ z + (1|group), 
                data_sample = samp2, 
                data_out_of_sample = non.samp2,
                method_y = "EBLUP",
                method_pi_score =  "EBLUP")
head(ht_EBLUP$tau_hat)

#My implementation leads to the same results as Setareh's implementation


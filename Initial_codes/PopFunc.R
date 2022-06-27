pop_data_generator <- function(ni_size=5,Ni_size=100,NoSim=1,
                               sigma2u=3, # variance of random effect outcome model 
                               sigma2e=6, #variance of the residuals outcome model
                               sigma2B = 1, # variance of the treatment effects
                               treat.mean=10,  # mean of treatment effects
                               sigma2V = 0.25, # variance of random effects for treatment status W
                               m=50, # number of areas
                               out.u=0.2, #percentace of outliers in u
                               out.mean.u= 9, # mean for outlier in u
                               out.v.u= 20, #variance for outliers in u
                               out.e=0.03,#percentace of outliers in e
                               out.mean.e= 20,
                               out.v.e= 150,
                               seed.start=2# seed to simulation with
                                )
{
  #ni <- rep(ni_size,m)
  Ni <- rep(Ni_size,m)
  N <- sum(Ni)
  #n <- sum(ni)
  pop.data <- vector("list",NoSim)
  for(j in 1:NoSim){
   set.seed(seed.start*j)
   u=rnorm(m,0,sqrt(sigma2u))
    
    ## outlier in area random effect
    k1 <- as.integer(out.u * m)
    u1 <- rnorm(k1, out.mean.u, sqrt(out.v.u))
    uu <- u
    uu[(m-k1+1):m] <- u1
    out.u.1 <- rep(out.u, m)
    u <- ifelse(out.u.1 > 0, uu, u)
    u=rep(u,times=Ni)
    
    ## Outlier in the residuals
    n1<-rbinom(N,1,out.e)
    mean.e<-20
    e <- (1-n1)*rnorm(N, 0, sqrt(sigma2e))+n1*rnorm(N, out.mean.e, sqrt(out.v.e))
    gr=rep(1:m,times=Ni)
    
    ## random effect for treatment assignment 
    V <- rnorm(m,0,sqrt(sigma2V))
    V <- rep(V, times=Ni)
    
    ## Treatment status
    X=matrix(c(rlnorm(N,log(4.5)-0.5,0.5)),nrow=N,ncol=1)
    Z=matrix(runif(N,min=0,max=1),nrow=N,ncol=1)
    pi<-exp(-1+0.5*Z+V)/(1+exp(-1+0.5*Z+V))
    W<-rbinom(N,1,pi)
    
    ## Treatment effect
    B <- rnorm(m, treat.mean, sd=sqrt(sigma2B))
    B <- rep(B, times=Ni)
    
    ##final outcome
    y=100+2*X+1*Z+ B*W +u+e # generating a scenario with different area specific effect
    
    pop.matrix<-cbind(y, X, Z, gr, W, pi)
    pop<-as.data.frame(pop.matrix)
    names(pop)<-c("y","x","z","area","w", "true_pi")    
    pop.data[[j]] <- pop 
  }
  
  return(pop.data)
}
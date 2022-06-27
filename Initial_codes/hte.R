##
#devtools::install_github("FelixSkarke/mquantreg")
#'
#'
#'
#'
#' @param formula_y
#' @param formula_treat
#' @param data Data table with 
#' @param method_outcome
#' @param method_pi_score
#'
#' @examples 
#' 
#' 
# Sketch
# Users provide population data + survey data
# We merge them 
populations <- generate_population(params = list(no_sim = 1, 
                                                 no_cov_unif = 5), 
                                   regression_type = "continuous",
                                   start_seed = 1)
samples <- generate_sample(populations, ni_size = 5)
#' tau_true <- calculate_tau(populations)
#'

hte <- function(formula_y, formula_treat, 
                data = sample_data, method_outcome,
                method_pi_score){
  
  outcome_regression
  
  pi_score_regression 
  
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

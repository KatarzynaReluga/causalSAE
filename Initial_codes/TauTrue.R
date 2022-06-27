tau_calc <- function(pop.data # The output of popFunc
                     ) {
  
  NoSim = length(pop.data)
  m <- length(unique(pop.data[[1]]$area))
  tau_true <- matrix(data = NA, nrow = NoSim , ncol = m)
  for(j in 1:NoSim){
    pop <- pop.data[[j]]
    ar=unique(pop$area)
    for(oo in ar){
    tau_true[j,oo] <- (sum(pop$w[pop$area==oo]*pop$y[pop$area==oo]/pop$true_pi[pop$area==oo])/sum(pop$w[pop$area==oo]/pop$true_pi[pop$area==oo]))-(sum((1-pop$w[pop$area==oo])*pop$y[pop$area==oo]/(1-pop$true_pi[pop$area==oo]))/sum((1-pop$w[pop$area==oo])/(1-pop$true_pi[pop$area==oo])))
    }
  }
  return(tau_true)
}
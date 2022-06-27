pop_sample <- function(pop.data,ni.size=5,samp=TRUE # to select the sampled part if true and non sampled part if not true 
                       ){
  
  NoSim = length(pop.data)
  m <- length(unique(pop.data[[1]]$area))
  samp.data <- vector("list",NoSim)
  non.samp.data <- vector("list", NoSim)
  ni <- rep(ni.size,m)

  for(j in 1:NoSim){
    pop <- pop.data[[j]]
    s <- strata(pop,"area", size=ni , method = "srswor") 
    samp.data[[j]] <- pop[s[,2],]
    non.samp.data[[j]] <-pop[-s[,2],]
  }
  if(samp){return(samp.data)}else{return(non.samp.data)}
}
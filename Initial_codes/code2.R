## Check if X is passed
ind.number = length(real.betas) - 1
if (any(X == FALSE)) {
  if(ind.number == 0){
    X = rep(1, n)
    dim(X) = c(n, 1)
    X = data.frame(X)
    colnames(X) = c("intercept")
  }
  else{
    X = c(rep(1, n), rnorm(n * ind.number, mean = 0, sd = 1))
    dim(X) = c(n, length(real.betas))
    X = data.frame(X)
    colnames(X) = c("intercept", paste0("x", 1:ind.number))
  }
  
}
data = as.data.frame(X)
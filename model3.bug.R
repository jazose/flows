model {
  for(i in 1:n){
    f[i] ~ dpois(lambda[i])
    lambda[i] <- exp(alpha1*log(hist[i]+1)+beta)
  }
  alpha1 ~ dunif(0,10)
  beta ~ dunif(-100,100)
}
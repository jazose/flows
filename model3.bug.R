model {
  for(i in 1:n){
    f[i] ~ dpois(lambda[i])
    lambda[i] <- exp(alpha1*log(hist[i]+1)+alpha2)
  }
  alpha1 ~ dunif(0,100)
  alpha2 ~ dunif(-100,100)
}
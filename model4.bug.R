model {
  for(i in 1:n){
    f[i] ~ dpois(lambda[i])
    lambda[i] <- exp(alpha1*log(hist[i]+1)+alpha2*log(o[i])+alpha3*log(d[i])+beta)
  }
  alpha1 ~ dunif(0,10)
  alpha2 ~ dunif(-10,10)
  alpha3 ~ dunif(-10,10)
  beta ~ dunif(-100,100)
}
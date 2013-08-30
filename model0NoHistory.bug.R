model {
  for(i in 1:n){
    f[i] ~ dpois(lambda[i])
    lambda[i] <- exp(alpha2*log(o[i])+alpha3*log(d[i])+
                       #delta1*dist1[i]+delta2*dist2[i]+delta3*dist3[i]+delta4*dist4[i]+
                       #delta5*dist5[i]+delta6*dist6[i]+
                       delta7*dist7[i]+
                       epsilon[i])
    epsilon[i] ~ dnorm(beta, tau)
  }
  alpha2 ~ dunif(-10,10)
  alpha3 ~ dunif(-10,10)
  #delta1 ~ dunif(-10,10)
  #delta2 ~ dunif(-10,10)
  #delta3 ~ dunif(-10,10)
  #delta4 ~ dunif(-10,10)
  #delta5 ~ dunif(-10,10)
  #delta6 ~ dunif(-10,10)
  delta7 ~ dunif(-10,10)
  beta ~ dunif(-10,10)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,10)
}
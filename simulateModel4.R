#Given input in the form of draws from the posterior distributions for Model 3,
#generate a sequence of simulated flows.

setup=function(){
  
  if(file.exists("C:/Users/jonazose")){
    setwd("C:/Users/jonazose/Dropbox/RA/Code/flows/flows_git")
  }else if(file.exists("C:/Users/Jon-Everyday")){
    setwd("C:/Users/Jon-Everyday/Dropbox/RA/Code/flows/flows_git/")
  }else if(file.exists("C:/Users/Jon")){
    setwd("C:/Users/Jon/Dropbox/RA/Code/flows/flows_git/")
  }else{
    setwd("/homes/jonazose/RA/flows_git/flows/")
  }
  
  #Read in data
  rawDat=scan("./Abel/flows_Decadal.txt")
  countryList<<-scan("./Abel/countryList_Decadal.txt",what="",sep=" ")
  flowArray<<-array(rawDat[4:length(rawDat)],dim=rawDat[1:3])
  rm(rawDat)
  rawPopDat=scan("./Abel/popDatMatrix_Decadal.txt");
  popDatMatrix<<-matrix(rawPopDat[3:length(rawPopDat)],nrow=rawPopDat[1],ncol=rawPopDat[2])
  
  #Toss out some countries that aren't in the CEPII database.
  #Guam (GUM), Mayotte (MYT), and US Virgin Islands (VIR)
  tossOutIndices=which(countryList %in% c("GUM","MYT","VIR"))
  flowArray<<-flowArray[-tossOutIndices,-tossOutIndices,];
  shortCountryList<<-countryList[-tossOutIndices];
  popDatMatrix<<-popDatMatrix[-tossOutIndices,]
  nm=length(shortCountryList);
  
  #####################
  #Read in CEPII stuff#
  #####################
  
  if(!file.exists("distanceArray.txt")){#If we didn't already make the distance array
    distData=read.csv("dist_cepii.csv",header=TRUE)
    modifiedCountryList=shortCountryList;
    modifiedCountryList[which(shortCountryList=="COD")]="ZAR"#COD in WorldBank is ZAR in CEPII
    modifiedCountryList[which(shortCountryList=="TLS")]="TMP"#TLS in WorldBank is TMP in CEPII
    modifiedCountryList[which(shortCountryList=="PSE")]="PAL"#PSE in WorldBank is PAL in CEPII
    
    distanceArray<<-array(0,dim=c(nm,nm,12))
    cat("Constructing distance matrix\n")
    for(i in 1:nm){
      cat(i,"\n")
      for(j in 1:nm){
        distanceArray[i,j,]=as.numeric(distData[which(distData$iso_o==modifiedCountryList[i]
                                                      & distData$iso_d==modifiedCountryList[j]),
                                                3:14])
      }
    }
    write(dim(distanceArray),"distanceArray.txt")
    write(distanceArray,"distanceArray.txt",append=TRUE)
  }else{#If we did make the distance array already, just read it in.
    distanceArrayDat=scan("distanceArray.txt")
    distanceArray<<-array(distanceArrayDat[4:length(distanceArrayDat)],
                          dim=distanceArrayDat[1:3]);
  }
  
  #Convert everything to vector form
  vectorLength=nm*(nm-1);#Keep track of the length of a single year's worth of data
  flowMatrix<<-matrix(0,nrow=dim(flowArray)[3],ncol=vectorLength);#Construct a matrix
  #where each row is a single year's data
  for(i in 1:dim(flowArray)[3]){
    M=flowArray[,,i];
    flowMatrix[i,]=M[row(M)!=col(M)]
  }
  flowMatrix<<-flowMatrix;
  
  originMatrix<<-matrix(rep(shortCountryList,nm),nrow=nm,byrow=FALSE);
  originVector<<-originMatrix[row(originMatrix)!=col(originMatrix)];
  destMatrix<<-t(originMatrix);
  destVector<<-destMatrix[row(destMatrix)!=col(destMatrix)];
  
  distanceMatrix<<-matrix(0,nrow=12,ncol=vectorLength);
  for(i in 1:12){
    M=distanceArray[,,i];
    distanceMatrix[i,]=M[row(M)!=col(M)];
  }
  distanceMatrix<<-distanceMatrix;
}

setup();

#Use for decadal data
y=as.vector(flowMatrix[c(2,3,4),]);
x=as.vector(flowMatrix[c(1,2,3),]);

#Construct vectors of origin and destination populations at beginning of the decade
oVec=numeric(0);
dVec=numeric(0);
for(j in 1:3){
  temp=rep(0,length(originVector));
  for(i in 1:length(originVector)){
    countryIndex=which(shortCountryList==originVector[i])
    temp[i]=popDatMatrix[countryIndex,j]
  }
  oVec=c(oVec,temp);
  
  temp=rep(0,length(destVector));
  for(i in 1:length(destVector)){
    countryIndex=which(shortCountryList==destVector[i])
    temp[i]=popDatMatrix[countryIndex,j]
  }
  dVec=c(dVec,temp);
}

generateFlows4=function(history, paramFile){
  #Generate predictions and return a matrix of medians and 80% and 95% confidence bounds
  
  #First, read in the posterior draws for alpha1 and alpha2.
  paramData=scan(paramFile);
  #Split into components
  nPosteriorDraws=length(paramData)/4;
  alpha1=paramData[1:nPosteriorDraws];
  alpha2=paramData[(nPosteriorDraws+1):(2*nPosteriorDraws)];
  alpha3=paramData[(2*nPosteriorDraws+1):(3*nPosteriorDraws)];
  beta=paramData[(3*nPosteriorDraws+1):(4*nPosteriorDraws)];
  
  #Initialize a matrix for simulation results
  simulatedFlows=matrix(0,nrow=nPosteriorDraws,ncol=length(history));
  for(i in 1:nPosteriorDraws){
    temp=rep(0,length(history));
    lambdaVec = exp(alpha1[i]*log(history+1)+alpha2[i]*log(oVec)+alpha3[i]*log(dVec)+beta[i])
    temp=rpois(n=length(history),lambda=lambdaVec)
    simulatedFlows[i,]=temp;
  }
  
  resultMatrix=matrix(0,nrow=length(x),ncol=5);
  for(i in 1:length(x)){
    resultMatrix[i,]=quantile(simulatedFlows[,i],c(0.025,0.1,0.5,0.9,0.975))
  }
  
  return(resultMatrix);
}

r=generateFlows4(x, "./Output/model4Output.txt")


#Return a mean absolute value of the median predictions minus the true values
MAE=mean(abs(y-r[,3]));
cat("MAE =",MAE);
MAE_log=mean(abs(log(y+1)-log(r[,3]+1)));
cat("MAE_log =",MAE_log)
maxAbsError=max(abs(y-r[,3]));
cat("Max error =",maxAbsError);

#Interval scores
#80% I.S.
is80=sum(r[,4]-r[,2]);
for(i in 1:nrow(r)){
  if(y[i]<r[i,2]){is80=is80+2/0.2*(r[i,2]-y[i]);}
  if(y[i]>r[i,4]){is80=is80+2/0.2*(y[i]-r[i,4]);}  
}
is80=is80/length(y)
cat("80% Interval score =",is80)

#95% I.S.
is95=sum(r[,5]-r[,1]);
for(i in 1:nrow(r)){
  if(y[i]<r[i,1]){is95=is95+2/0.2*(r[i,1]-y[i]);}
  if(y[i]>r[i,5]){is95=is95+2/0.2*(y[i]-r[i,5]);}  
}
is95=is95/length(y)
cat("95% Interval score =",is95)

#######################
#Confidence interval inclusion
#######################
#80% CI
inclusion80=(y>=r[,2])&(y<=r[,4])
cat("80% Inclusion =",mean(inclusion80))
inclusion95=(y>=r[,1])&(y<=r[,5])
cat("95% Inclusion =",mean(inclusion95))
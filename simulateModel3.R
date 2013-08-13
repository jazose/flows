#Given input in the form of draws from the posterior distributions for Model 3,
#generate a sequence of simulated flows.

setup=function(){
  
  if(file.exists("C:/Users/jonazose")){
    setwd("C:/Users/jonazose/Dropbox/RA/Code/flows/flows_git")
  }else if(file.exists("C:/Users/Jon-Everyday")){
    setwd("C:/Users/Jon-Everyday/Dropbox/RA/Code/flows/flows_git/")
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

generateFlows3=function(history, paramFile){
  #Generate predictions and return a matrix of medians and 80% and 95% confidence bounds
  
  #First, read in the posterior draws for alpha1 and alpha2.
  paramData=scan(paramFile);
  #Split into components
  nPosteriorDraws=length(paramData)/2;
  alpha1=paramData[1:nPosteriorDraws];
  beta=paramData[(nPosteriorDraws+1):(2*nPosteriorDraws)];
  
  #Initialize a matrix for simulation results
  simulatedFlows=matrix(0,nrow=nPosteriorDraws,ncol=length(history));
  for(i in 1:nPosteriorDraws){
    temp=rep(0,length(history));
    lambdaVec = exp(alpha1[i]*log(history+1)+beta[i])
    temp=rpois(n=length(history),lambda=lambdaVec)
    simulatedFlows[i,]=temp;
  }
  
  return(simulatedFlows);
}

simulatedFlows=generateFlows3(x, "./Output/model3Output.txt")
medianVec=rep(0,length(x));
for(i in 1:length(x)){
  medianVec[i]=median(simulatedFlows[,i])
}

#Return a mean absolute value of the median predictions minus the true values
MAE=mean(abs(y-medianVec));
cat("MAE =",MAE);
MAE_log=mean(abs(log(y+1)-log(medianVec+1)));
cat("MAE_log =",MAE_log)
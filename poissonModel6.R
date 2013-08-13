#Run a simulated Poisson hierarchical model.
#Model 5 takes log(lambda) to be a linear combination of log(previous flow + 1), log(pop. at origin), and log(pop. at destination),
#   and also the various distance metrics.
#Just like model 5, but includes a couple more distance metrics.

library(rjags);library(coda);

###############
#Read in data.#
###############

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

#Use for annual data
#y=as.vector(flowMatrix[c(15,25,35),])
#x=as.vector(flowMatrix[c(5,15,25),])

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


#Choose a small sample from the y and x vectors
sampleFraction=1 #Use all of the data for fitting
compressedDataSize=floor(sampleFraction*length(y))
compressedDataIndices=sort(sample.int(n=length(y),size=compressedDataSize))
smallY=y[compressedDataIndices];
smallX=x[compressedDataIndices];

##################
#Run through JAGS#
##################

#The compressed version
jags <- jags.model('model6.bug.R',
                   data=list('n' = compressedDataSize,
                             'hist' = smallX,
                             'f' = smallY,
                             'o' = oVec,
                             'd' = dVec,
                             'dist1' = rep(distanceMatrix[1,],3),
                             'dist2' = rep(distanceMatrix[2,],3),
                             'dist3' = rep(distanceMatrix[3,],3),
                             'dist4' = rep(distanceMatrix[4,],3),
                             'dist5' = rep(distanceMatrix[5,],3),
                             'dist6' = rep(distanceMatrix[6,],3),
                             'dist7' = rep(distanceMatrix[7,],3)),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 750)

samples=coda.samples(jags,
                     c('alpha1','alpha2','alpha3','beta','delta1','delta2','delta3','delta4'),
                     1000,
                     thin=5)

sampleData=as.matrix(rbind(samples[[1]],samples[[2]],samples[[3]],samples[[4]]))
write(sampleData,"./Output/model6output.txt")
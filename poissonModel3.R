#Run a simulated Poisson hierarchical model.

library(rjags);library(coda);

if(file.exists("C:/Users/jonazose")){
  setwd("C:/Users/jonazose/Dropbox/RA/Code/flows/")
}else{
  setwd("C:/Users/Jon-Everyday/Dropbox/RA/Code/flows/")
}

###############
#Read in data.#
###############

setup=function(){
  
  if(file.exists("C:/Users/jonazose")){
    setwd("C:/Users/jonazose/Dropbox/RA/Code/flows/")
  }else{
    setwd("C:/Users/Jon-Everyday/Dropbox/RA/Code/flows/")
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


#Choose a small sample from the y and x vectors
sampleFraction=0.01 #Use only 1% of the data for fitting
compressedDataSize=floor(sampleFraction*length(y))
compressedDataIndices=sort(sample.int(n=length(y),size=compressedDataSize))
smallY=y[compressedDataIndices];
smallX=x[compressedDataIndices];

##################
#Run through JAGS#
##################

# #The uncompressed version
# jags <- jags.model('model2.bug.R',
#                    data = list('n' = length(y),
#                                'hist' = x,
#                                'f' = y),
#                    n.chains = 4,
#                    n.adapt = 5000)

#The compressed version
jags <- jags.model('model2.bug.R',
                   data=list('n' = compressedDataSize,
                             'hist' = smallX,
                             'f' = smallY),
                   n.chains = 4,
                   n.adapt = 1000)

update(jags, 10000)

samples=coda.samples(jags,
                     c('alpha1','alpha2'),
                     5000,
                     thin=1)

sampleData=as.matrix(rbind(samples[[1]],samples[[2]],samples[[3]],samples[[4]]))
plot(sampleData[,1])
plot(sampleData[,2])
hist(sampleData)
#write(sampleData,"model2Output_August6.txt")
#August 6 output used 1000 for adaptation, 200000 for update, and 5000 draws, thinned by 250.
plot(density(sampleData[,1]));mean(sampleData[,1]);
plot(sampleData[,1],sampleData[,2],xlab="alpha1",ylab="alpha2")





#####################
#Look at geometry of posterior density surface
#####################

#UNFINISHED

logLik=function(xVec,yVec,alpha1,alpha2){
  l=0;
  for(i in 1:length(xVec)){
    lambda=exp(alpha1*log(hist[i]+1)+alpha2*eps[i])
    l=l+
  }
  return(l);
}

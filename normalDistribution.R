
#packages need
require(mvtnorm)
require(randtoolbox)
require(numDeriv)

#specify model parameters
true=list()
true$mu=c(3,5)
true$sigma=matrix(c(1.0,0.9,
                    0.9,2.0),nrow=2,ncol=2)
cov2cor(true$sigma)

parNames=c("mu1","mu2")

#specify prior
prior=list()
prior$mu=c(1,1)
prior$sigma=matrix(c(10.0, 0.0,
                     0.0 , 10.0),nrow=2,ncol=2)

#generate data
nData=100
dat=rmvnorm(nData,true$mu,true$sigma)

par(mfrow=c(1,2))
hist(dat[,1],main="x1",col=rgb(1,0,0,.5),breaks=20,prob=T,xlim=range(as.numeric(dat)),ylim=c(0,max(c(density(dat[,1])$y,density(dat[,2])$y))+.1))
hist(dat[,2],main="x2",col=rgb(0,0,1,.5),breaks=20,prob=T,add=T)
plot(dat,xlab="x1",ylab="x2",pch=16,col=rgb(0,1,0,.5),main="generated data")

##################################
####### optimization parameters
#################################3
nParsModel=length(parNames)
nParsPost=2*nParsModel
S=10 #number of samples to approximate the KL divergence/ELBO
nChains=nParsPost*3 #number of particles for optimization
maxIter=600
noise=0.001 #uniform noise parameter for solution proposal
purifyProb=.20 
mutationProb=1


#################################
####### advanced parameters
#################################
negInf=-750
# qFunction=c("MFVB","LRVB")
# qFunction=qFunction[2]

useQMC=F
quasiSeq=c("sobol","halton")
quasiSeq=quasiSeq[1]

useStopCrit=F
stopCrit=-Inf
##################################



#####################################################
######## Model Functions 
#####################################################
Transform.Parameters=function(x){
  #pass labeled parameter vector and transforms the vector
  x[names(x)=="mu1"]<-x[names(x)=="mu1"]
  x[names(x)=="mu2"]<-x[names(x)=="mu2"]
  #x[names(x)=="sigma"]<-exp(x[names(x)=="sigma"])
  x=x[!is.na(x)]
  return(x)
}


Log.Prior=function(useTheta){
  require(mvtnorm)
  #returns log of prior density of theta sample
  out=0 #initalize output variable
  #assign parameter names
  names(useTheta)<-parNames
  useTheta<-Transform.Parameters(useTheta)
  
  out=dmvnorm(useTheta,prior$mu,prior$sigma,log=T)
  #cactch any NA or infinite values
  out[!is.finite(out)]<-negInf
  return(out)
}

Log.Like=function(useTheta,useData){
  #returns log likelihood of theta sample
  
  out=0 #initalize output variable
  
  #assign parameter names
  names(useTheta)<-parNames
  
  #transform parameters
  useTheta=Transform.Parameters(useTheta)
  
  #calculate log density
  logDens=dmvnorm(useData,mean=c(useTheta["mu1"],useTheta["mu2"]),true$sigma,log=T)
  
  #cactch any NA or infinite values
  logDens[!is.finite(logDens)]<-negInf
  
  #sum log densities
  out=sum(logDens)
  
  return(out)
  
}


#####################################################
######## VB Functions 
#####################################################

q.Log=function(useTheta,useLambda){
  #returns log density q(theta|lambda)
  out=0
  out=dnorm(x=useTheta[1:nParsModel],
            mean=useLambda[1:nParsModel],
            sd=exp(useLambda[(nParsModel+1):(nParsPost)]),
            log=T)
  out[(out==-Inf) | is.na(out)]=negInf
  return(sum(out));
}


q.Sample=function(useLambda,S){
  #returns a S by nParameter matrix sampled from q(theta|lambda)
  out=matrix(NA,S,nParsModel)
  if(useQMC==T){
    require(randtoolbox)
    if(quasiSeq=="sobol")quantileMat=sobol(S,nParsModel)
    if(quasiSeq=="halton")quantileMat=halton(S,nParsModel)
  } else {
    quantileMat=matrix(runif(S,0,1),S,nParsModel)
  }
  for(i in 1:nParsModel){
    qs=quantileMat[,i]
    out[,i]=qnorm(qs,useLambda[i],exp(useLambda[nParsModel+i]))
  }
  return(out)
}

KL.hat=function(lambda,data=dat){
  #monte carlo approximation KL divergence up to a constant
  
  out=0 #initalize output vector
  
  #sample from q
  thetaMat<-q.Sample(useLambda=lambda,S) 
  
  #calc mean differences in log densities for thetaMat
  for(s in 1:S){ #can parallelize/vector this 
    out=out+q.Log(useTheta = thetaMat[s,],useLambda = lambda)
    out=out-Log.Like(useTheta = thetaMat[s,],data) - Log.Prior(useTheta = thetaMat[s,]) 
  }
  
  out=out/S
  
  return(out)
}

####################################
#### run VB!
################################

#declare arrays
lambdaArr=array(NA,dim=c(nParsPost,nChains,maxIter))
weightArr=array(Inf,dim=c(nChains,maxIter))

#initalize arrays
xInit=rep(0,nParsPost)

for(i in 1:nChains){
  while(!is.finite(weightArr[i,1])){
    
    #sample values of lambda 
    lambdaArr[,i,1]=rnorm(nParsPost,xInit,.1) 
    
    #calculate weight
    weightArr[i,1]=KL.hat(lambda=lambdaArr[,i,1],dat) 
  }
  print(paste("initialized:", i,"/",nChains))
}

stopCritMet=FALSE
t=1

while( t<maxIter & (!stopCritMet) ){
  t=t+1
  
  for(i in 1:nChains){ #this is the part easily parallelized
    
    current <- lambdaArr[,i,t-1] 
    currentWeight <- weightArr[i,t-1]
    
    #######################################
    #purification
    #######################################
    if(runif(1)<purifyProb){
      temp <- KL.hat(lambda=current,dat)
      if(is.finite(temp)){
        currentWeight=temp
      }
    }
    
    #######################################
    #mutation & crossover
    #######################################
    if(runif(1)<mutationProb){
      
      #scaling factor
      gamma=runif(1,.5,1.5) 
      
      #sample 2 parent vectors to mate
      parentIdx=sample((1:nChains)[-i],size = 2,replace = F) 
      
      parent1=lambdaArr[,parentIdx[1],t-1]
      parent2=lambdaArr[,parentIdx[2],t-1]
      
      #cross parents and generate a proposal
      proposal = current + gamma*(parent1-parent2)+runif(1,-noise,noise)
      
      proposalWeight=KL.hat(lambda=proposal,dat)
      
      #accept or reject porposal  
      if(!is.finite(proposalWeight))proposalWeight=Inf
      if(proposalWeight<currentWeight){
        current <- proposal
        currentWeight <- proposalWeight
      }
      
    }
    #update array
    current -> lambdaArr[,i,t] 
    currentWeight -> weightArr[i,t]
    
  }
  if(t%%10==0){
    print(paste("iteration: ",t,"/",maxIter,sep=""))
  }
  
}
lastIter=t

#########################################
### plot the particle trajectories
#########################################
#graphical parameters

lty0=1
lwd0=1.5

par(mfrow=c(1,2))
matplot(apply(weightArr*-1,MARGIN = 2,median),type="l",main="ELBO median",xlab="iteration",ylab="ELBO median",lty=lty0,lwd=lwd0)
matplot(log(apply(weightArr*-1,MARGIN = 2,var)),type="l",main="ELBO variance ",xlab="iteration",ylab="log(ELBO variance)",lty=lty0,lwd=lwd0)


#########################################
### plot the particle trajectories
#########################################
#graphical parameters
lty0=1
lwd0=1.25

par(mfrow=c(2,2),mar=c(4,4,4,4))
parNamesTrans=parNames#c("\mu1","log(\sigma)")

for(k in 1:nParsPost){
  if(k>nParsModel){
    label=paste(parNamesTrans[k-nParsModel],"posteriorSD")
    matplot(t(exp(lambdaArr[k,,])),type="l",main=label,lty=lty0,lwd=lwd0,ylab="",xlab="iteration")
  } else {
    label=paste(parNamesTrans[k],"posteriorMean")
    matplot(t(lambdaArr[k,,]),type="l",main=label,lty=lty0,lwd=lwd0,ylab="",xlab="")
  }
}

#########################################
#### plot the posteriors
#########################################
meanPost=numeric(nParsModel)
covPost=diag(nParsModel)
nSamplesPost=1000

#########################################
#### samples from VB posts
lambdaOptimal=lambdaArr[,which.min(weightArr[,lastIter]),lastIter]
meanPost<-lambdaOptimal[1:nParsModel]
diag(covPost)<-exp(2*lambdaOptimal[(nParsModel+1):nParsPost])
postSamplesVB=rmvnorm(nSamplesPost,meanPost,covPost)

#########################################
#### samples from true post
postSigma=solve(solve(prior$sigma)+nData*solve(true$sigma))
postMean=postSigma %*% (solve(prior$sigma) %*% prior$mu + nData*solve(true$sigma)%*%apply(dat,MARGIN=2,FUN=mean))
postSamplesTrue=rmvnorm(nSamplesPost,postMean,postSigma)

#########################################
par(mfrow=c(2,2))

hist(postSamplesTrue[,1],main="mu1",col=rgb(0,0,1,.5),breaks=20,prob=T,xlim=range(as.numeric(c(postSamplesTrue[,1],postSamplesVB[,1]))),ylim=c(0,max(c(density(postSamplesVB[,1])$y,density(postSamplesTrue[,1])$y))+.1))
hist(postSamplesVB[,1],main="mu1",col=rgb(1,0,0,.5),breaks=20,prob=T,add=T)

hist(postSamplesTrue[,2],main="mu2",col=rgb(0,0,1,.5),breaks=20,prob=T,xlim=range(as.numeric(c(postSamplesTrue[,2],postSamplesVB[,2]))),ylim=c(0,max(c(density(postSamplesVB[,2])$y,density(postSamplesTrue[,2])$y))+.1))
hist(postSamplesVB[,2],main="mu2",col=rgb(1,0,0,.5),breaks=20,prob=T,add=T)

plot(postSamplesTrue,col=rgb(0,0,1,.15),xlab=parNames[1],ylab=parNames[2],pch=16,xlim = range(c(postSamplesVB[,1],postSamplesTrue[,1])),ylim = range(c(postSamplesVB[,2],postSamplesTrue[,2])))
points(postSamplesVB,col=rgb(1,0,0,.15),pch=16)

plot(NA,xlim=c(0,1),ylim=c(0,1),main = "legend",frame.plot=F,xlab="",ylab="",yaxt='n' ,xaxt='n', ann=FALSE)
legend("center", 
       c("true","VB"), fill=c("blue","red"), horiz=TRUE, cex=0.8)


#########################################
#### correct the posteriors!
#########################################
require(numDeriv)
useQMC=T
S=25
hessKL=hessian(KL.hat,lambdaOptimal)
if(!all(eigen(hessKL)$values > (-10e-6))) print("WARNING: optimization procedure failed")

LRVBcovMat=solve(hessKL)[1:nParsModel,1:nParsModel]
postSamplesVB=rmvnorm(nSamplesPost,meanPost,LRVBcovMat)

#########################################
par(mfrow=c(2,2))

hist(postSamplesTrue[,1],main="mu1",col=rgb(0,0,1,.5),breaks=20,prob=T,xlim=range(as.numeric(c(postSamplesTrue[,1],postSamplesVB[,1]))),ylim=c(0,max(c(density(postSamplesVB[,1])$y,density(postSamplesTrue[,1])$y))+.1))
hist(postSamplesVB[,1],main="mu1",col=rgb(1,0,0,.5),breaks=20,prob=T,add=T)

hist(postSamplesTrue[,2],main="mu2",col=rgb(0,0,1,.5),breaks=20,prob=T,xlim=range(as.numeric(c(postSamplesTrue[,2],postSamplesVB[,2]))),ylim=c(0,max(c(density(postSamplesVB[,2])$y,density(postSamplesTrue[,2])$y))+.1))
hist(postSamplesVB[,2],main="mu2",col=rgb(1,0,0,.5),breaks=20,prob=T,add=T)

plot(postSamplesTrue,col=rgb(0,0,1,.15),xlab=parNames[1],ylab=parNames[2],pch=16,xlim = range(c(postSamplesVB[,1],postSamplesTrue[,1])),ylim = range(c(postSamplesVB[,2],postSamplesTrue[,2])))
points(postSamplesVB,col=rgb(1,0,0,.15),pch=16)

plot(NA,xlim=c(0,1),ylim=c(0,1),main = "legend",frame.plot=F,xlab="",ylab="",yaxt='n' ,xaxt='n', ann=FALSE)
legend("center", 
       c("true","VB"), fill=c("blue","red"), horiz=TRUE, cex=0.8)


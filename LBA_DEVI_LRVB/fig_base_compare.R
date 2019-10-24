#load("demcmc.RData")
#mainDir<-setwd("~/Dropbox/18_projects/LikelihoodFreeVarBayes/ProofOfConcept/GD vs DE")
start=10
#################################colors here
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

alpha=1
alpha=.1
alpha_hist=.75
mc_col=add.alpha("#2c3e50",1)
MFVB_col=add.alpha("#3498db",alpha)
LRVB_col=add.alpha("#e74c3c",alpha)


mc_colh=add.alpha("#ecf0f1",1)
MFVB_colh=add.alpha("#3498db",alpha_hist)
LRVB_colh=add.alpha("#e74c3c",alpha_hist)
################################################


setwd(file.path(mainDir, subDir))
tnmc=length(keep.samples)

if(plot.lower==TRUE)theta=array(NA,c(n.chains,n.pars,tnmc))
if(plot.weights==TRUE)weights=array(NA,c(n.chains,tnmc))

for(q in 1:n.chains){
  if(plot.lower==TRUE){
    
    temp=t(as.matrix(read.table(paste("chain",q,"_lower.txt",sep=""),header=F)))
    theta[q,,]=temp[,1:tnmc]
  }
  
  if(plot.weights==TRUE){
    temp=t(as.matrix(read.table(paste("chain",q,"_weights.txt",sep=""),header=F)))
    weights[q,]=temp[,1:tnmc]
  }
  print(round(q/n.chains*100))
}

#######################################################################################

breaks=50

if(plot.lower==TRUE){
  par(mfrow=c(2,2),ask=F)
  for(k in 1:n.pars){
    #if(any(k==c(1,2,3,6)))abline(h=log(true$vec[k]),lwd=3)
    #if(any(k==c(4,5)))abline(h=logit(true$vec[k]),lwd=3)
    if(par.names[k]=="v1"){
      matplot(t(exp(theta[,k,seq(start,tnmc,by=extra.thin)])),type="l",lty=1,ylab=k)
      
      hist(exp(theta[,k,seq(start,tnmc,by=extra.thin)]),prob=T,breaks=breaks,xlab=k, main = par.names[k])
      abline(v=(true$v[1]),lwd=2)
    }  
    if(par.names[k]=="v2"){
      matplot(t(exp(theta[,k,seq(start,tnmc,by=extra.thin)])),type="l",lty=1,ylab=k)
      hist(exp(theta[,k,seq(start,tnmc,by=extra.thin)]),prob=T,breaks=breaks,xlab=k, main = par.names[k])
      abline(v=true$v[2],lwd=2)
    }  
    
    if(par.names[k]=="b"){
      matplot(t(exp(theta[,k,seq(start,tnmc,by=extra.thin)])),type="l",lty=1,ylab=k)
      hist(exp(theta[,k,seq(start,tnmc,by=extra.thin)]),prob=T,breaks=breaks,xlab=k, main = par.names[k])
      abline(v=true$b,lwd=2)
    }  
    
    if(par.names[k]=="z"){
      matplot(t(invlogit(theta[,k,seq(start,tnmc,by=extra.thin)])),type="l",lty=1,ylab=k)
      hist(invlogit(theta[,k,seq(start,tnmc,by=extra.thin)]),prob=T,breaks=breaks,xlab=k, main = par.names[k])
      abline(v=true$z,lwd=2)
    }  
    
    if(par.names[k]=="tau"){
      matplot(t(invlogit(theta[,k,seq(start,tnmc,by=extra.thin)])),type="l",lty=1,ylab=k)
      hist(invlogit(theta[,k,seq(start,tnmc,by=extra.thin)]),prob=T,breaks=breaks,xlab=k, main = par.names[k])
      abline(v=true$tau,lwd=2)
    }  
    
    #if(any(k==c(1,2,3,6)))abline(v=log(true$vec[k]),lwd=3,col="red")
    #if(any(k==c(4,5)))abline(v=logit(true$vec[k]),lwd=3,col="red")
  }}


if(plot.weights==TRUE){
  par(mfrow=c(1,1),ask=T)
  
  matplot(t(weights[,seq(start,tnmc,by=extra.thin)]),type="l",lty=1,main=paste("Subject"),ylab="log likelihood")
}



#############################
tr.theta=theta
for(k in 1:n.pars){
  
  
  
  if(par.names[k]=="tau" | par.names[k]=="z"){tr=invlogit(theta[,k,seq(start,tnmc,by=extra.thin)])
  }else{tr=exp(theta[,k,seq(start,tnmc,by=extra.thin)])
  }
  
  tr.theta[,k,seq(start,tnmc,by=extra.thin)]=tr  
  
}



#######################################
par(mfrow=c(5,5),ask=F,mar=c(1,1,1,1,1))
theta_seq=seq(start,tnmc,by=extra.thin)
size=150
ss=sample(theta_seq,size)
size2=900

####################################
require(mvtnorm)

#extract lambda parameters
muMFVB=numeric(n.pars)
for(k in 1:n.pars){
  muMFVB[k]=mean(lambda.arr[(last-10):last,,k])
}
sigmaMFVB=numeric(n.pars)
for(k in 1:n.pars){
  sigmaMFVB[k]=exp(2*mean(lambda.arr[(last-10):last,,k+n.pars]))
}

#sample MFVB
MFVBsamples<-rmvnorm(size*n.chains,muMFVB, diag(sigmaMFVB))
MFVBsamples2<-rmvnorm(size2*n.chains,muMFVB, diag(sigmaMFVB))


#transform MFVB samples
for(k in 1:n.pars){
  if(par.names[k]=="tau" | par.names[k]=="z"){
    MFVBsamples[,k]=invlogit(MFVBsamples[,k])
      MFVBsamples2[,k]=invlogit(MFVBsamples2[,k])
  }else{MFVBsamples[,k]=exp(MFVBsamples[,k])
  MFVBsamples2[,k]=exp(MFVBsamples2[,k])
  }
}

LRVBsamples<-rmvnorm(size*n.chains,muMFVB, LRCovMat)

LRVBsamples2<-rmvnorm(size2*n.chains,muMFVB, LRCovMat)
#transform MFVB samples
for(k in 1:n.pars){
  if(par.names[k]=="tau" | par.names[k]=="z"){LRVBsamples[,k]=invlogit(LRVBsamples[,k])
  LRVBsamples2[,k]=invlogit(LRVBsamples2[,k])
  }else{LRVBsamples[,k]=exp(LRVBsamples[,k])
  LRVBsamples2[,k]=exp(LRVBsamples2[,k])
  }
}
require(MASS)
require(akima)
ymax=c(14,12,25,30,25)
xmin=c(.95,.5,.75,.3,.375)
xmax=c(1.5,1.05,1.15,.65,.69)
breaks=c(120,80,100,100,70)
par(mfrow=c(2,3),ask=F)

######################################### plot it
for(i in 1:k){
  for(j in 1:k){
    #plot(tr.theta[,i,ss],tr.theta[,j,ss],xlab=par.names[i],ylab=par.names[j],col=mc_col)    
    
    if(i==j){
    hist(as.numeric(tr.theta[,i,]),prob=T,main="",ylim=c(0,ymax[i]),xlim=c(xmin[i],xmax[i]),col=mc_colh,breaks=breaks[i],xlab="density")
    lines(density(LRVBsamples2[,i]),col=LRVB_colh,lwd=2.5)
    lines(density(MFVBsamples2[,i]),col=MFVB_colh,lwd=2.5)}
  #   }else if(i>j){plot(NA,xlim=c(0,1),ylim=c(0,1),yaxt="n",xaxt="n",frame.plot=FALSE)
  #     text(.5,.4,cex=1.25,col=LRVB_colh,paste(as.character(round(cov2cor(LRCovMat)[i,j],4))))
  #     text(.5,.6,cex=1.25,paste(as.character(round(cor(as.numeric(theta[,i,]),as.numeric(theta[,j,])),4))))
  #   }else{
  #   f<-kde2d(as.numeric(tr.theta[,j,ss]),as.numeric(tr.theta[,i,ss]),n=20,h=.12)
  #   contour(f, drawlabels = F,col=mc_col)
  #   points(LRVBsamples[,j],LRVBsamples[,i],col=LRVB_col,pch=16)
  #   points(MFVBsamples[,j],MFVBsamples[,i],col=MFVB_col,pch=16)
  #   } 
  # }
}
}


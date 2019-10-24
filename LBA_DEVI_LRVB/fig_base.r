
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


# 
# count=1
# if(plot.hyper==TRUE){
# par(mfrow=c(2,2),ask=T)
# for(k in 1:n.hpars){
# if(plot.priors==TRUE){
# if(k==(n.pars+1))count=1
# xs=seq(min(phi[,k,start:tnmc]),max(phi[,k,start:tnmc]),length=200)
# if(k<=n.pars)ys=dnorm(xs,prior[[count]]$mu,prior[[count]]$sigma)
# if(k>n.pars)ys=dinvgamma(xs,prior[[count]]$alpha,prior[[count]]$beta)
# }
# matplot(t(phi[,k,start:tnmc]),type="l",lty=1,main="",ylab=hpar.names[k])
# hist(phi[,k,start:tnmc],prob=T,breaks=breaks,main="",xlab=hpar.names[k])
# if(plot.priors==TRUE){
# lines(xs,ys,lty=2)
# count=count+1
# }}}

if(plot.weights==TRUE){
par(mfrow=c(1,1),ask=T)

matplot(t(weights[,seq(start,tnmc,by=extra.thin)]),type="l",lty=1,main=paste("Subject"),ylab="log likelihood")
}

setwd(file.path(mainDir))


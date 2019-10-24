#########################################################################
log.dens.prior=function(theta){
  dens=0
  require(msm)
  names(theta)<-par.names
  dens=dens+sum(dnorm(as.numeric(theta["v1"]),0,.5,log=T))
  dens=dens+sum(dnorm(as.numeric(theta["v2"]),0,.5,log=T))
  dens=dens+sum(dnorm(as.numeric(theta["b"]),0,.5,log=T))
  dens=dens+sum(dnorm(as.numeric((theta["z"])),0,1/4,log=T))
  dens=dens+sum(dnorm(as.numeric((theta["tau"])),0,1/4,log=T))
  return(dens)
}

#specify log likelihood function
log.like.dens=log.dens.like=function(use.theta,use.data){
  ##### VARIABLE DESCRIPTION
  out=0
  v0=c(exp(use.theta[1]),exp(use.theta[2]))
  b0=exp(use.theta[3])
  z0=invlogit(use.theta[4])
  tau0=min(use.data$rt)*invlogit(use.theta[5])
  
  like=get.dens.2choice(rt=use.data$rt,resp=use.data$resp,b=b0,A=b0*z0,v=v0,s=true$s,t0=tau0)
  like=log(like)
  like[(like==-Inf) | is.na(like)]=-Inf
  out=sum(like)
  
  
  return(sum(out));
}

require(compiler)
log.like.dens<-cmpfun(log.like.dens)
log.dens.like=log.like.dens



log.like=function(use.theta,use.data){
  ##### VARIABLE DESCRIPTION
  like=0
  v0=c(exp(use.theta[1]),exp(use.theta[2]))
  b0=exp(use.theta[3])
  z0=invlogit(use.theta[4])
  tau0=min(use.data$rt)*invlogit(use.theta[5])
  
  like=get.dens.2choice(rt=use.data$rt,resp=use.data$resp,b=b0,A=b0*z0,v=v0,s=true$s,t0=tau0)
  like=log(like)
  like[(like==-Inf) | is.na(like)]=-Inf
  like=sum(like)
  
  return((like));
}


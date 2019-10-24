crossover=function(i,pars,use.theta,use.like,use.data){
  use.weight=use.like[i]
  gamma = runif(1,.5,1)#2.38/sqrt(2*length(pars))

  index=sample(c(1:n.chains)[-i],2,replace=F)
  theta=use.theta[i,]						
  theta[pars]=use.theta[i,pars] + gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) + runif(1,-b,b)
  theta=matrix(theta,1,length(theta))
  like=log.dens.like(theta,use.data)
  weight=like + log.dens.prior(theta)
  if(is.na(weight))weight=-Inf
  if(runif(1) < exp(weight-use.weight)) {							
    use.theta[i,]=theta
    use.like[i]=like
  }
  c(use.like[i],use.theta[i,])
}


migrate=function(pars,use.theta,use.like,use.data){
  lnum1=sample(c(1:n.chains),1)						# determine how many groups to work with
  lnum2=sort(sample(c(1:n.chains),lnum1,replace=F))			# which groups we will work with
  thetaset=matrix(NA,lnum1,length(use.theta[1,]))					# initialize
  currentset=propset=propw=currw=numeric(lnum1)
  index=numeric(lnum1)
  for(i in 1:lnum1){
    index[i]=sample(1:n.chains,1,replace=F)	
    thetaset[i,]=use.theta[lnum2[i],]
    thetaset[i,pars]=use.theta[lnum2[i],pars] + runif(1,-b,b)			# create a set of these particles to swap
    currentset[i]=use.like[lnum2[i]]
    currw[i]=currentset[i] + log.dens.prior(use.theta[lnum2[i],])
  }
  if(lnum1==1)propset=log.dens.like(thetaset[1,],data)
  if(lnum1!=1)propset=unlist(sfLapply(1:lnum1,function(x,data,theta)log.dens.like(theta[x,],data),theta=thetaset,data=data))
  propw=propset + apply(thetaset,1,log.dens.prior)
  if(runif(1) < exp(propw[lnum1] - currw[1])){
    use.theta[lnum2[1],]=thetaset[lnum1,]					# swap the first with the last (creating a circle)
    use.like[lnum2[1]]=propset[lnum1]
  }
  if(lnum1!=1){								# make sure we are not done yet
    for(i in 1:(lnum1-1)){		
      if(runif(1) < exp(propw[i] - currw[i+1])){
        use.theta[lnum2[i+1],]=thetaset[i,]					# swap the first with the last (creating a circle)
        use.like[lnum2[i+1]]=propset[i]
      }}}
  list(weight=use.like,theta=use.theta)
}
purify_chains=function(i,pars,use.theta,use.like,use.data){
  theta=use.theta[i,]
  like=log.dens.like(theta,use.data)
  weight=like + log.dens.prior(theta)
  if(!is.na(weight) & weight!=-Inf){
    use.like[i]=like    
  }
  c(use.like[i],use.theta[i,])
}


write.files=function(q,use.theta,use.weight,append=TRUE){
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
  }
  write(round(use.theta[q,],6),paste("chain",q,"_lower.txt",sep=""),ncolumns=n.pars,append=append)
  write(round(use.weight[q],8),paste("chain",q,"_weights.txt",sep=""),ncolumns=1,append=append)
  setwd(mainDir)
}

########################################## run it
sfInit(parallel=TRUE, cpus=cores, type="SOCK")
sfClusterSetupRNG()
theta=array(NA,c(n.chains,n.pars))
weight=array(-Inf,c(n.chains))
sfExportAll(except=list("theta","weight"))

for(i in 1:n.chains){

    while(weight[i]==-Inf){
      theta[i,]=rnorm(n.pars,x.init,.5)
      weight[i]=log.dens.like(theta[i,],use.data=data) + log.dens.prior(theta[i,])
      if(is.na(weight[i]))weight[i]=-Inf
    }
  print(paste("Initialization ",round(i/n.chains*100),"% Complete",sep=""))
}

junk=sfLapply(1:n.chains,write.files,use.theta=theta,use.weight=weight,append=FALSE)

########################################## 

for(i in 2:nmc){

  temp=matrix(unlist(sfLapply(1:n.chains,crossover,pars=1:n.pars,use.theta=array(theta,c(n.chains,n.pars)),use.like=weight,use.data=data)),n.chains,n.pars+1,byrow=T)
  weight=temp[,1]
  theta=temp[,2:(n.pars+1)]
   
  if(i<migrate.duration){
    if(runif(1)<migrate.prob){
      out=migrate(use.theta=theta,use.like=weight,use.data=data)
      weight=out$weight
      theta=out$theta
    }}
  if(i%%purify.rate==0){
    temp=matrix(unlist(sfLapply(1:n.chains,purify_chains,pars=1:n.pars,use.theta=array(theta,c(n.chains,n.pars)),use.like=weight,use.data=data)),n.chains,n.pars+1,byrow=T)
    weight=temp[,1]
    theta=temp[,2:(n.pars+1)]
  }
  
  if(any(i==keep.samples))temp=sfLapply(1:n.chains,write.files,use.theta=theta,use.weight=weight,append=T)
  if(i%%10==0)print(i)
}

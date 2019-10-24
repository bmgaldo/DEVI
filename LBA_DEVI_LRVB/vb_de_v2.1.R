#########################################################################
#vectorized q.loq
require("compiler")
logit=function(x){
  log(x/(1-x))
}

invlogit=function(x){
  1/(1+exp(-x))
}

#DE Functions
require(snowfall)
require(MASS)
require(msm)

crossover_vb=function(i,pars,use.lambda,use.kl,data,time){
  #proposes new parameters values for chain i using a linear function of the other chains
  #################################################################
  #pars refer to which posterior paramaters you're updating
  #use.lambda is the current posterior parameters 
  #use.kl is the current fitness (the negative KL divergence)
  #data is a list of data
  #################################################################################
  use.weight=use.kl[i]
  gamma = 2.38/sqrt(length(pars))
  index=sample(c(1:n.chains)[-i],2,replace=F)
  lambda=use.lambda[i,]						
  lambda[pars]=use.lambda[i,pars] + gamma*(use.lambda[index[1],pars]-use.lambda[index[2],pars]) + runif(1,-b,b)
  lambda=matrix(lambda,1,length(lambda))
  kl=ELBO.hat(lambda,data,S=S,t=time)
  weight=kl 
  if(is.na(weight))weight=-Inf
  if(weight>use.weight){							
    use.lambda[i,]=lambda
    use.kl[i]=kl
  }
  c(use.kl[i],use.lambda[i,])
}


purify_chains_vb=function(i,pars,use.lambda,use.kl,data,time){
  #re-evaluates the fitness of a chain
  #prevents chain sticking since the objective function is stochastic.
  lambda=use.lambda[i,]
  kl=ELBO.hat(lambda,data,S,t=time)
  weight=kl
  if(!is.na(weight) & weight!=-Inf){
    use.kl[i]=kl    
  }
  c(use.kl[i],use.lambda[i,])
}
migrate_vb=function(use.lambda,use.kl,log.dens,S,data){
  pars=dim(use.lambda)[2]
  lnum1=sample(c(1:n.chains),1)										# determine how many groups to work with
  lnum2=sort(sample(c(1:n.chains),lnum1,replace=F))							# which groups we will work with
  lambdaset=matrix(NA,lnum1,pars)									# initialize
  currentset=propset=propw=currw=numeric(lnum1)
  index=numeric(lnum1)
  for(i in 1:lnum1){
    index[i]=sample(1:n.chains,1,replace=F)	
    lambdaset[i,]=use.lambda[lnum2[i],] + runif(1,-b,b)				# create a set of these particles to swap
    propset[i]=log.dens(lambdaset[i,],data,S)
    currentset[i]=use.kl[lnum2[i]]
    propw[i]=propset[i]
    currw[i]=currentset[i]
  }
  if(propw[lnum1] > currw[1]){
    use.lambda[lnum2[1],]=lambdaset[lnum1,]							# swap the first with the last (creating a circle)
    use.kl[lnum2[1]]=propset[lnum1]
  }
  if(lnum1!=1){											# make sure we are not done yet
    for(i in 1:(lnum1-1)){		
      if(propw[i] > currw[i+1]){
        use.lambda[lnum2[i+1],]=lambdaset[i,]							# swap the first with the last (creating a circle)
        use.kl[lnum2[i+1]]=propset[i]
      }}}
  list(weight=use.kl,lambda=use.lambda)
}


#########################################################################
#VB Functions
#specify Q, paramteric distribution used to approximate posterior

q.log=function(theta,use.lambda,j){
  out=NULL
  out=dnorm(theta[1:n.pars.model],use.lambda[1:n.pars.model],exp(use.lambda[(n.pars.model+1):n.pars.post]),log=T)
  out[(out==-Inf) | is.na(out)]=neg.inf
  return(sum(out));
}

q.log.deriv=function(theta,use.lambda,j){
  out=NULL
  out=dnorm(theta[1:n.pars.model],use.lambda[1:n.pars.model],(use.lambda[(n.pars.model+1):n.pars.post]),log=T)
  out[(out==-Inf) | is.na(out)]=neg.inf
  return(sum(out));
}



#sample from q so u can estimate Kullback-Leibler divergence at these points
q.sample=function(use.lambda,S,t){
  require("msm")
  require(mvtnorm)
  require(randtoolbox)
  out=matrix(NA,S,n.pars.model)
  #qs.fix=sobol(S,5)
  
  for(i in 1:n.pars.model){
    qs=runif(S,0,1)#qs.fix[,i]+runif(S,-.001,.001)
    temp=qnorm(qs,use.lambda[i],exp(use.lambda[n.pars.model+i]))
    out[,i]=temp
  }
  return(out);
}

#samples q for KL derivative approximation
q.sample.deriv=function(use.lambda,S,t){
  require("msm")
  require(mvtnorm)
  require(randtoolbox)
  out.mat=matrix(NA,S,n.pars.model)
  qs.fix=sobol(S,5)
  
  for(i in 1:n.pars.model){
    qs=qs.fix[,i]
    temp=qnorm(qs,use.lambda[i],(use.lambda[n.pars.model+i]))
    out.mat[,i]=temp
  }
  return(out.mat);
}

#approximates log negative KL divergence 
ELBO.hat=function(use.lambda,data=data,S=S,t=1){
  
  theta.mat<-q.sample(use.lambda,S,t)
  out=0
  for(i in 1:S){
    out=out+q.log(theta.mat[i,],use.lambda)-(log.like(use.theta=theta.mat[i,],use.data=data)+log.dens.prior(theta=theta.mat[i,]))
  }
  out=(out*-1)/S
  
  return(out)
}


#approximates KL for numerical derivative approximators
KLforDeriv=function(use.lambda,data=data.list,S=10,t=1){
  
  theta.mat<-q.sample.deriv(use.lambda,S,t)
  out=0
  for(i in 1:S){
    out=out+q.log.deriv(theta.mat[i,],use.lambda)-(log.like(use.theta=theta.mat[i,],use.data=data)+log.dens.prior(theta=theta.mat[i,]))
  }
  out=(out)/S
  
  return(out)
}



########################################## run it
sfInit(parallel=TRUE, cpus=cores, type="SOCK")
sfClusterSetupRNG()


##################################intialize
##################################
  lambda.arr=array(NA,c(iters,n.chains,n.pars.post))
  weight.arr=array(-Inf,c(iters,n.chains))
  for(i in 1:n.chains){
    while(weight.arr[1,i]==-Inf){
      lambda.arr[1,i,]=rnorm(n.pars.post,x.init.vb,.5)
      weight.arr[1,i]=ELBO.hat(use.lambda=lambda.arr[1,i,],data=data,S)
    }
    print(paste(round(i/n.chains*100),"% initialization"))
  }

  sfExportAll(except=list("lambda.arr","weight.arr"))

  lb.hat=Inf
  i=2
  

while(abs(lb.hat>tol) & (i<=iters)){
  temp=matrix(unlist(sfLapply(1:n.chains,crossover_vb,pars=1:n.pars.post,time=i,use.lambda=array(lambda.arr[i-1,,],c(n.chains,n.pars.post)),use.kl=weight.arr[i-1,],data=data)),n.chains,n.pars.post+1,byrow=T)
  weight.arr[i,]=temp[,1]
  lambda.arr[i,,]=temp[,2:(n.pars.post+1)]
  ##################################
  if(i%%purifyTimeVB==0){
    temp=matrix(unlist(sfLapply(1:n.chains,purify_chains_vb,pars=1:n.pars.post,time=i,use.lambda=array(lambda.arr[i,,],c(n.chains,n.pars.post)),use.kl=weight.arr[i,],data=data)),n.chains,n.pars.post+1,byrow=T)
    weight.arr[i,]=temp[,1]
    lambda.arr[i,,]=temp[,2:(n.pars.post+1)]
  }
  
  ###################################
  if(runif(1)<migrate.prob){
    temp=migrate_vb(use.lambda=array(lambda.arr[i,,],c(n.chains,n.pars.post)),use.kl=weight.arr[i,],log.dens=ELBO.hat,S=S,data=data)
    weight.arr[i,]=temp$weight
    lambda.arr[i,,]=temp$lambda
  }
  ##################################
  if(i%%10==0)print(i)

  i=i+1
}


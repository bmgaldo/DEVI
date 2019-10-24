###############
###############
###############
###############
rm(list=ls())
setwd("~/Dropbox/18_projects/LikelihoodFreeVarBayes/ProofOfConcept/LBA_DEVA_LRVB copy")
mainDir=getwd()
########################################## load the functions you will use

require(snowfall)
require(MASS)
require(msm)

source("lba.r") #model functions

###########################################################
logit=function(x){
  log(x/(1-x))
}

invlogit=function(x){
  1/(1+exp(-x))
}

#######################################
##### Generate Data to fit
######################################
n.choices=2
n=750
true=NULL #create list for true parameters
true$b=1 #Decision Boundary
true$z=.5 #bias parameter
true$A=true$b*true$z #start point
true$v=c(1.2,.8) #drift rate
true$s=1 #between trial variability
true$t0=.15 #nondecision
data=rlba(n=n,b=true$b,A=true$A,vs=true$v,s=true$s,t0=true$t0)
shift=0
data$log.rt=log(data$rt+shift)
true$tau=true$t0/min(data$rt)


#initialization vectors
par.names=c("v1","v2","b","z","tau")
n.pars=length(par.names)
x.init=c(log(true$v),log(true$b),logit(true$z),logit(true$tau))

#load the priors and likelihood function
source("inference_functions.R")
################################################################################

ptm <- proc.time()

######################################################################
#######################################
##### MCMC - AlGORITHM PARAMETERS
######################################
n.chains=30
n.chains.demcmc=n.chains
nmc=5500
burnin=500
thin=4
b=.001 #uniform noise parameter, crossover step
migrate.duration=round(.75*burnin)
migrate.prob=.05
purify.rate=5
keep.samples=seq(burnin,nmc,thin)
cores=4
N=10000
neg.inf=-750

subDir="MCMC"
#########################################
source("de_mcmc.r")

time.demcmc <- proc.time() - ptm
save.image("demcmc.RData")

setwd(mainDir)
start=10

extra.thin=1
start.weights=10
setwd(mainDir)
plot.lower=T
plot.weights=T

source("fig_base.R")

######################################################################
setwd(mainDir)
ptm <- proc.time()
#######################################
##### VB AlGORITHM PARAMETERS
######################################
likelihoodFree=F
iters=500#when to give up, avoid infinite loop and set finite memory
n.pars.model=length(x.init) #number of parameters in the generative model
n.pars.post=2*n.pars.model #number of parameters in the function approximating the posterior
S=6
#number of samples to approximate KL divergence
n.chains=30#n.pars.post*3
purify.rate=5
shuffle.prob=0
migrate.prob=0.1 #
purifyTimeVB=5
b=.001 #uniform noise parameter, crossover step
tol=-Inf#10^-5
Rvalue=0
M=1
neg.inf=-750
N=round(10000)
firstRun=T
x.init.vb=c(x.init,rep(log(.25),n.pars))
#########################################
source("vb_de_v2.1.r")
last=i-1
time.devb = proc.time() - ptm
saveRDS(lambda.arr,"lambdaVB")
saveRDS(weight.arr,"weightVB")
saveRDS(last,"lastVB")
saveRDS(time.devb,"timeVB")


###########################################plotting
lambda.arr<-readRDS("lambdaVB")
weight.arr<-readRDS("weightVB")

iters=500
#get estimates
wind=10
lambda.est=numeric(n.pars.post)
for(i in 1:n.pars.post){
  lambda.est[i]=mean(lambda.arr[(iters-wind):iters,,i])
}

start=1
tnmc=iters
extra.thin=1
par(mfrow=c(2,2), ask=F)
for(i in 1:(n.pars.post)){
  if(i>n.pars.model){
    matplot(((lambda.arr[,,i])),type="l",lty=1, main=paste("posterior log(sd)",par.names[i-n.pars]) )
    #    abline(h=log((true$post.var[i-n.pars.model])^.5),col="black",lwd=2)
  }
  if(i<=n.pars.model){
    matplot((lambda.arr[,,i]),type="l",lty=1,main=paste("posterior mean",par.names[i]))
    #abline(h=true$post.mean[i],col="black",lwd=2)
  }
  
  abline(h=lambda.est[i],lwd=2)
}
matplot((weight.arr)/n,main="weight",type="l",lty=1)


############################################ correct pars
############################################
require(numDeriv)
plot.lower=T
plot.weights=T

data.list<-data
x=lambda.est
x[6:10]<-exp(x[6:10])
hess=hessian(KLforDeriv,x)
print(round(eigen(hess)$values,4))
hess.inv=solve(hess)
#plot(diag(hess.inv^.5)[1:5],exp(lambda.est[6:10]))
newSDs=diag(hess.inv^.5)[1:n.pars.model]
LRCovMat<-hess.inv[1:n.pars.model,1:n.pars.model]

setwd(mainDir)
source("fig_base_compare.r")


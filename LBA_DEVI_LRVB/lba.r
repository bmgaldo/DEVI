


fptcdf=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max==0) return(pnorm(chi/z,mean=driftrate,sd=sddrift,lower.tail=F))
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
  chizu=chiminuszu/zs ; chizumax=xx/zs
  tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
  tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
  1+(tmp1+tmp2)/x0max
}

fptpdf=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max==0) return( (chi/z^2)*dnorm(chi/z,mean=driftrate,sd=sddrift))
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
  (driftrate*(pnorm(chizu)-pnorm(chizumax)) + 
     sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max
}


rlba=function(n,b,A,vs,s,t0){
  rt=resp=numeric(n)
  choices=length(vs)
  for(i in 1:n){
    d=0
    while(all(d<=0)){
      d=rnorm(choices,vs,s)
    }
    k=runif(choices,0,A)
    ttf=(b-k)/d
    min.rt=min(ttf[ttf>0])
    rt[i]=min.rt + t0
    resp[i]=c(1:choices)[ttf==min.rt]
  }
  list(rt=rt,resp=resp)
}

######################################

raw=function(x,choices,b,A,vs,s,t0){
  d=0
  while(any(d<=0)){
    d=rnorm(choices,vs,s)
  }
  k=runif(choices,0,A)
  ttf=(b-k)/d
  ttf + t0
}

rlba_raw=function(n,b,A,vs,s,t0){
  resp=numeric(n)
  choices=length(vs)
  t(sapply(1:n,raw,choices,b,A,vs,s,t0))
}


#######################################


n1PDF=function(t,A,b,v,sdv) {
  N=length(v) # Number of responses.
  if (N>2) {
    tmp=array(dim=c(length(t),N-1))
    for (i in 2:N) tmp[,i-1]=fptcdf(t=t,A=A,b=b,v=v[i],sdv=sdv)
    G=apply(1-tmp,1,prod)    
  } else {
    G=1-fptcdf(t=t,A=A,b=b,v=v[2],sdv=sdv)
  }
  G*fptpdf(t=t,A=A,b=b,v=v[1],sdv=sdv)
}

n1PDFfixedt0=function(t,x0max,chi,drift,sdI,truncdrifts=TRUE) {
  # Generates defective PDF for responses on node #1.
  # "truncdrifts" sets whether that part of the multi-variate
  # normal distribution on drift rates which would otherwise
  # lead to non-terminating trials is truncated.
  N=length(drift) # Number of responses.
  if (N>2) {
    tmp=array(dim=c(length(t),N-1))
    for (i in 2:N) tmp[,i-1]=fptcdf(z=t,x0max=x0max[i],chi=chi[i],
                                    driftrate=drift[i],sddrift=sdI[i])
    G=apply(1-tmp,1,prod)
  } else {
    G=1-fptcdf(z=t,x0max=x0max[2],chi=chi[2],driftrate=drift[2],sddrift=sdI[2])
  }
  out=G*fptpdf(z=t,x0max=x0max[1],chi=chi[1],driftrate=drift[1],sddrift=sdI[1])
  if (truncdrifts) {
    out=out/(1-prod(pnorm(-drift/sdI)))
    out[t<=0]=0
    return(out)
  } else {
    return(out)
  }
}

get.dens.2choice=function(rt,resp,b,A,v,s,t0){
  sapply(1:length(rt),function(x,rt,resp,b,A,s,v,t0)n1PDFfixedt0(rt[x]-t0,x0max=A,chi=b,drift=c(v[resp[x]],v[-resp[x]]),sdI=s),b=c(b,b),A=c(A,A),s=c(s,s),v=v,t0=t0,resp=resp,rt=rt)
}

get.dens.2choice=function(rt,resp,b,A,v,s,t0){
  out=numeric(length(rt))
  b=as.numeric(b)
  A=as.numeric(A)
  v=as.numeric(v)
  t0=as.numeric(t0)
  tmp=(resp==1)
  out[tmp] =n1PDFfixedt0(rt[tmp]-t0,x0max=c(A,A),chi=c(b,b),drift=v,s=c(s,s))
  out[!tmp]=n1PDFfixedt0(rt[!tmp]-t0,x0max=c(A,A),chi=c(b,b),drift=v[2:1],s=c(s,s))
  out
}



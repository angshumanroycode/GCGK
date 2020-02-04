gcgk.preset<-function(d,parameters=NULL,B.parameters=1000,normalize="normalized"){
  if(is.null(parameters)){
    D=dist(replicate(d,runif(B.parameters)))
    shl=quantile(D,c(0.5,0.05))
    m=ceiling(log2(shl[1]/shl[2]))
    parameters=shl[1]/2^(0.5+0:m)
  }
  L=list()
  L[[1]]=parameters
  pl=length(parameters)
  if("normalized" %in% normalize){
    C=numeric(pl)
    for(i in 1:pl){
      s=parameters[i]
      v1=V(0,1,0,1,s/sqrt(d))
      l<-function(x) (sqrt(2*pi)*s*(pnorm(x/s)+pnorm((1-x)/s)-1))^d
      v2=integrate(l,0,1)$value
      v3=(V(0,1,0,1,s))^d
      C[i]=1/sqrt(v1-2*v2+v3)
    }
    L[[2]]=C
  }
  return(L)
}

gcgk<-function(X,parameters=NULL,B.parameters=1000,
               normalize=c("unnormalized","normalized"),preset=NULL){
  if(is.null(preset))
    preset=gcgk.preset(nrow(X),parameters,B.parameters,normalize)
  parameters=preset[[1]]
  tstat=GCGK(X,parameters)
  result=list()
  resultnames=NULL
  l=1
  result[[l]]=parameters
  resultnames[l]="parameters"
  if("unnormalized" %in% normalize){
    l=l+1
    result[[l]]=tstat
    resultnames[l]="unnormalized"
  }
  if("normalized" %in% normalize){
    l=l+1
    result[[l]]=preset[[2]]*tstat
    resultnames[l]="normalized"
  }
  names(result)=resultnames
  return(result)
}

gcgk.test<-function(X,method=c("sum","max","fdr","nsum","nmax"),
                    B=1000,parameters=NULL,B.parameters=1000,
                    alpha=0.05,preset=NULL){
  if("nsum" %in% method || "nmax" %in% method){
    normalize="normalized"
  } else {
    normalize="unnormalized"
  }
  if(is.null(preset))
    preset=gcgk.preset(nrow(X),parameters,B.parameters,normalize)
  parameters=preset[[1]]
  L=GCGKTEST(X,parameters,B)
  tstat=L[[1]]
  ndist=L[[2]]
  result=list()
  resultnames=NULL
  l=1
  if("sum" %in% method){
    tstat.sum=sum(tstat)
    ndist.sum=rowSums(ndist)
    pvalue=(sum(tstat.sum<ndist.sum)+1)/(length(ndist.sum)+1)
    reject=pvalue<alpha
    subresult=list(reject,pvalue)
    names(subresult)=c("reject","p.value")
    result[[l]]=subresult
    resultnames[l]="sum"
    l=l+1
  }
  if("max" %in% method){
    tstat.max=max(tstat)
    ndist.max=apply(ndist,1,max)
    pvalue=(sum(tstat.max<ndist.max)+1)/(length(ndist.max)+1)
    reject=pvalue<alpha
    subresult=list(reject,pvalue)
    names(subresult)=c("reject","p.value")
    result[[l]]=subresult
    resultnames[l]="max"
    l=l+1
  }
  if("fdr" %in% method){
    np=NCOL(ndist)
    nr=NROW(ndist)
    tstat.fdr=matrix(rep(tstat,nr),nrow=nr,byrow=T)
    pvalsort=sort(colMeans(tstat.fdr<ndist))
    Levels=alpha*(1:np)/np
    reject=sum(pvalsort<Levels)>0
    subresult=list(reject)
    names(subresult)=c("reject")
    result[[l]]=subresult
    resultnames[l]="fdr"
    l=l+1
  }
  if("normalized" %in% normalize){
    C=preset[[2]]
    tstat=t(t(L[[1]])*C)
    ndist=t(t(L[[2]])*C)
    if("nsum" %in% method){
      tstat.sum=sum(tstat)
      ndist.sum=rowSums(ndist)
      pvalue=(sum(tstat.sum<ndist.sum)+1)/(length(ndist.sum)+1)
      reject=pvalue<alpha
      subresult=list(reject,pvalue)
      names(subresult)=c("reject","p.value")
      result[[l]]=subresult
      resultnames[l]="nsum"
      l=l+1
    }
    if("nmax" %in% method){
      tstat.max=max(tstat)
      ndist.max=apply(ndist,1,max)
      pvalue=(sum(tstat.max<ndist.max)+1)/(length(ndist.max)+1)
      reject=pvalue<alpha
      subresult=list(reject,pvalue)
      names(subresult)=c("reject","p.value")
      result[[l]]=subresult
      resultnames[l]="nmax"
    }
  }
  names(result)=resultnames
  return(result)
}
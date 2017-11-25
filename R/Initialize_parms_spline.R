Ispline=function(intknots,lowept,upperept){
  k=3
  #determine knot sequence
  knotseq=c(rep(lowept,k+1),intknots,rep(upperept,k+1))
  numbases=length(knotseq)-k-2
  
  #create matrix of bases
  bases=matrix(NA,nrow=numbases,ncol=2)
  for(i in 1:numbases) bases[i,]=c(knotseq[i+1],knotseq[i+k+1])
  
  return(list(bases=bases,knotseq=knotseq))
}


getIsplineC=function(xtrue,knotseq,bases){
  
  numBases=nrow(bases)
  lxtrue=length(xtrue)
  
  mat=rep(0,numBases*lxtrue)
  
  storage.mode(mat) <- "double"
  storage.mode(knotseq) <- "double"
  storage.mode(xtrue) <- "double"
  storage.mode(numBases) <- "integer"
  storage.mode(lxtrue) <- "integer"
  temp=.C("getIspline",xtrue,lxtrue,knotseq,mat,numBases)
  designMatrix=matrix(temp[[4]],ncol=numBases,nrow=lxtrue)
  return(designMatrix)
  
}

initialize_parms_spline=function(xobs,yobs,xcensl,xcensu,xgrid){
  
  ## Initialize xtrue  
  xtrue=rep(NA,length(xobs))
  
  n=sum(xcensl==0 & xcensu==0)
  n1=sum(xcensl==1)
  n2=sum(xcensu==1)
  
  xtrue[xcensl==0 & xcensu==0]=xobs[xcensl==0 & xcensu==0]-runif(n,0,1)
  xtrue[xcensl==1]=xobs[xcensl==1]-runif(n1,.5,2)
  xtrue[xcensu==1]=xobs[xcensu==1]+runif(n2,.5,1.5)
  
  ### endpoints
  xlower=min(xobs)-2
  if(sum(xcensl)>0) xlower=min(xobs)-4
  xupper=max(xobs)+1
  if(sum(xcensu)>0) xupper=max(xobs)+3
  
  ### Initalize Isplines
  lowept=min(xobs)-.5
  upperept=max(xobs)+.5
  if(sum(xcensu==1)>0)
    upperept=max(xobs)+3
  if(sum(xcensl==1)>0)
    lowept=min(xobs)-4
  dist=1
  parms=Ispline(seq(min(xobs),max(xobs),by=dist),lowept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq

  #### Initial spline coefficients
  designMatrix=getIsplineC(xtrue,knotseq,bases)
  designMatrixGrid=getIsplineC(xgrid,knotseq,bases)
  
  min_sav=1e7
  for(i in seq(.5,7,by=.1)){
    icoefs1=seq(.1,i,length=ncol(designMatrix))
    ytrue1=as.numeric(icoefs1%*%t(designMatrix))
    if(sum(abs(yobs-ytrue1))<min_sav){ min_sav=sum(abs(yobs-ytrue1)); save=i}
  }
  coefs=seq(.1,save,length=ncol(designMatrix))

  return(list(coefs=coefs,xtrue=xtrue,designMatrix=designMatrix,designMatrixGrid=designMatrixGrid,
              xlower=xlower,xupper=xupper))
}

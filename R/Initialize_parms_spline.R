  initialize_parms_spline=function(xobs,yobs,xcensl,xcensu,xgrid){
  
  ## Initial xtrue
  xtrue=rep(NA,length(xobs))
  for(i in 1:length(xobs)){
    if(xcensl[i]==0 & xcensu[i]==0){
      xtrue[i]=xobs[i]-runif(1,0,1)
    }else if(xcensl[i]==1 & xcensu[i]==0){
      xtrue[i]=xobs[i]-runif(1,.5,1.5)
    }else if(xcensl[i]==0 & xcensu[i]==1){
      xtrue[i]=xobs[i]+runif(1,0,1)
    }
  }
  
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

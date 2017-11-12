initialize_parms_log=function(xobs,yobs,xcensl,xcensu){
  
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
  
  #### Initial Logistic coefficients
  nls_logistic=function(xtrue,coef){
    mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
    fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
    ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
        (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
       exp(coef[4]*(coef[2]-xtrue)))
    return(ytrue)
  }
  coefs=tryCatch({coef(nls(yobs~nls_logistic(xtrue,coef),start=list(coef=c(max(yobs),.2,.1,.1)),
          lower=c(0.01,0.01,0.01,0.01),algorithm='port',trace=FALSE,
          control=list(maxiter=1000,tol = 1e-04)))},
          error=function(cond){return(c(max(yobs),.2,.1,.1))})

  return(list(coefs=coefs,xtrue=xtrue))
}

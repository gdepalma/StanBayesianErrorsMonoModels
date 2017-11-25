initialize_parms_log=function(xobs,yobs,xcensl,xcensu){
  
  ## Initialize xtrue  
  xtrue=rep(NA,length(xobs))

  n=sum(xcensl==0 & xcensu==0)
  n1=sum(xcensl==1)
  n2=sum(xcensu==1)
  
  xtrue[xcensl==0 & xcensu==0]=xobs[xcensl==0 & xcensu==0]-runif(n,0,1)
  xtrue[xcensl==1]=xobs[xcensl==1]-runif(n1,.5,2)
  xtrue[xcensu==1]=xobs[xcensu==1]+runif(n2,.5,1.5)
  
  #### Initial Logistic coefficients
  nls_logistic=function(coef1,coef2,coef3,coef4){
    mb = (2*coef3*coef4)/(coef3+coef4)
    fx = 1/(1+exp(-mb*(coef2-xtrue)))
    ytrue=coef1*(fx*exp(coef3*(coef2-xtrue))+(1-fx)*exp(coef4*
        (coef2-xtrue)))/(1+fx*exp(coef3*(coef2-xtrue))+(1-fx)*
                             exp(coef4*(coef2-xtrue)))
    return(ytrue)
  }

  fit=tryCatch({
    nlsLM(yobs~nls_logistic(coef1,coef2,coef3,coef4),start=list(coef1=max(yobs)-2,coef2=.2,coef3=.1,coef4=.1))
  },
    error=function(cond) {
      message(cond)
      stop("Logistic initial parameters not converging.")
    }
  )
  
  return(list(coefs=coef(fit),xtrue=xtrue))
}

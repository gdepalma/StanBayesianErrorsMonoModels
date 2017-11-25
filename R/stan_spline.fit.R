
stan_spline.fit <- function(dat_sav,xgrid,nchains=1,xsig=.707, ysig=2.121) {
  
  xobs=dat_sav$xobs
  yobs=dat_sav$yobs
  xcensl=dat_sav$xcensl
  xcensu=dat_sav$xcensu
  ycensl=dat_sav$ycensl
  ycensu=dat_sav$ycensu
  Ngrid=length(xgrid)
  

  ### Initialize
  parms=initialize_parms_spline(xobs,yobs,xcensl,xcensu,xgrid)
  coefs=parms$coefs
  designMatrix=parms$designMatrix
  designMatrixGrid=parms$designMatrixGrid
  xtrue=parms$xtrue
  xlower=parms$xlower
  xupper=parms$xupper
  
  ### Run Stan Model
  rstan_options(auto_write = TRUE)
  options(mc.cores=parallel::detectCores())
  
  
  dat=list(y=yobs,xgrid=xgrid,Ngrid=Ngrid,n_groups=5,N=N,xobs=xobs,xsig=xsig, ysig=ysig,
           xcensl=xcensl,xcensu=xcensu,ycensl=ycensl,ycensu=ycensu,numCoef=length(coefs),
           designMatrix=designMatrix,designMatrixGrid,xlower=xlower,xupper=xupper)
  init_fun <- function() {list(mu=seq(min(xtrue)+1,max(xtrue)-1,length=5),sigma=rep(1,5),
            Theta=rep(1/5,5),xtrue=xtrue,coef=coefs,alpha_coef=1)}
  stanfit <- stanmodels$spline
  fit <- rstan::sampling(stanfit, pars = c("MIC_Dens","gx"),data = dat, iter = 500,chains = nchains,thin=5,init=init_fun)
  parms <- extract(fit)
  
  return(parms)
}

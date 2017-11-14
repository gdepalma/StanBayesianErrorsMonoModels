

#' @rdname stan_logistic
#' @export

stan_logistic.fit <- function(dat_sav,xgridk,nchains=1) {
  
  xobs=dat_sav$xobs
  yobs=dat_sav$yobs
  xcensl=dat_sav$xcensl
  xcensu=dat_sav$xcensu
  ycensl=dat_sav$ycensl
  ycensu=dat_sav$ycensu
  

  ### Initialize
  parms=initialize_parms_log(xobs,yobs,xcensu,xcensl)
  coefs=parms$coefs
  xtrue=parms$xtrue
  
  ### Run Stan Model
  options(mc.cores=parallel::detectCores())
  
  dat=list(y=yobs,xgrid=xgrid,Ngrid=Ngrid,n_groups=5,ysig=ysig,N=N,xobs=xobs,xsig=xsig,
           xcensl=xcensl,xcensu=xcensu,ycensl=ycensl,ycensu=ycensu)
  init_fun <- function() {list(mu=seq(min(xtrue)+1,max(xtrue)-1,length=5),sigma=rep(1,5),
            Theta=rep(1/5,5),xtrue=xtrue,coef=coefs)}
  stanfit <- stanmodels$logistic
  fit <- rstan::sampling(stanfit, pars = c("MIC_Dens","gx"),data = dat, iter = 500,chains = nchains,thin=2,init=init_fun)
  parms <- rstan::extract(fit)
  
  return(parms)
}

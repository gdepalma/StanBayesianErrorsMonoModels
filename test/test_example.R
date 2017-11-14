library(mixtools)
library(ggthemes)
library(grid)
library(gridExtra)
library(StanBayesianErrorsMonoModels)

### Generate Sim Data for four-parameter logistic model
nobs=400
xcens=rep(0,nobs)
ycens=rep(0,nobs)

xsig=.707; ysig=2.121
xgrid=seq(-12,12,length=1000)

popmn=c(-4.6,-2,1); popstd=c(1.1,1.5,1.5); popprob=c(.6,.2,.2)
xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
xobs=ceiling(xtrue+rnorm(nobs,0,xsig))

coef=c(35,1.17,.1,1.2)
mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
  (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
  exp(coef[4]*(coef[2]-xtrue)))
yobs=round(ytrue+rnorm(nobs,0,ysig))



### Set up data
xgrid=seq(min(xobs)-1,max(xobs)+1,length=1000)
Ngrid=length(xgrid)
N=length(xobs)
xcensl=rep(0,N)
xcensl[xcens==-1] = 1
xcensu=rep(0,N)
xcensu[xcens==1] = 1
ycensl=rep(0,N)
ycensl[ycens==-1] = 1
ycensu=rep(0,N)
ycensu[ycens==1] = 1

dat_sav=data.frame(xobs,yobs,xcensl,xcensu,ycensu,ycensl)

list_of_draws = stan_logistic.fit(dat_sav,xgrid,nchains=1)

parms=output_graphs(list_of_draws,xgrid,dat_sav)
pltRel=parms$pltRel
pltDens=parms$pltDens

plt1 <- ggplot_gtable(ggplot_build(pltRel))
plt2 <- ggplot_gtable(ggplot_build(pltDens))
maxWidth = unit.pmax(plt1$widths[2:3], plt2$widths[2:3])
plt1$widths[2:3] <- maxWidth
plt2$widths[2:3] <- maxWidth
plot(grid.arrange(plt1, plt2, ncol=1, heights=c(5,2)))


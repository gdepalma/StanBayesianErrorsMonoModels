
library(StanBayesianErrorsMonoModels)

a1 = read_csv(file='test/data1.csv')
xobs=a1$MIC; yobs=a1$DIA; xcens=rep(0,length(xobs)); ycens=rep(0,length(yobs))

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


# list_of_draws = stan_logistic.fit(dat_sav,xgrid,nchains=1)

list_of_draws = stan_spline.fit(dat_sav,xgrid,nchains=1,numIter=500)

list_of_draws = stan_spline.fit(dat_sav,xgrid,nchains=3,numIter=500)


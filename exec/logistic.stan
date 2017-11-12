functions{
	real getytrue(vector coef, real xtrue,real yobs){
	
		real mb;
		real fx;
		real ytrue;
	
		mb = (2*coef[3]*coef[4])/(coef[3]+coef[4]);
		fx = 1/(1+exp(-mb*(coef[2]-xtrue)));
		ytrue = coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
			(coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
			exp(coef[4]*(coef[2]-xtrue)));
			
		return(ytrue);

	}
	
  real xlikecens(real xobs, real xtrue,real xsig,real xcensu,real xcensl){

	  if(xcensu==0 && xcensl==0)
	     return(log_diff_exp(normal_lcdf(xobs | xtrue,xsig),normal_lcdf(xobs-1 | xtrue,xsig)));
	  else if(xcensl==0 && xcensu==1){
	     return(log_diff_exp(0,normal_lcdf(xobs-1|xtrue,xsig)));
	  }
	  else if(xcensl==1 && xcensu==0)
	     return(normal_lcdf(xobs | xtrue,xsig));
	  else return(0);
	}
	
	real ylikecens(real yobs, real ytrue,real ysig,real ycensu,real ycensl){
  
	  if(ycensu==0 && ycensl==0)
	     return(log_diff_exp(normal_lcdf(yobs+.5 | ytrue,ysig),normal_lcdf(yobs-.5 | ytrue,ysig)));
	  else if(ycensl==0 && ycensu==1)
	     return(log_diff_exp(0,normal_lcdf(yobs-.5|ytrue,ysig)));
	  else if(ycensl==1 && ycensu==0)
	     return(normal_lcdf(yobs+.5 | ytrue,ysig));
	  else return(0);
	}
}

data {
  int N;
  vector[N] yobs;
  vector[N] xobs;
  vector[N] xcensl;
  vector[N] xcensu;
  vector[N] ycensl;
  vector[N] ycensu;
  int n_groups;
  int Ngrid;
  real xgrid[Ngrid];
  real xsig;
  real ysig;
}

parameters {
  real mu[n_groups]; 
  simplex[n_groups] v;
  real<lower=0> sigma[n_groups]; 
  vector<lower=0>[4] coef;
  real xtrue[N];
}

transformed parameters{
  simplex [n_groups] Theta;
  Theta[1] = v[1];
  // stick-break process based on The BUGS book Chapter 11 (p.294)
  for(j in 2:(n_groups-1)){
      Theta[j]= v[j]*(1-v[j-1])*Theta[j-1]/v[j-1]; 
  }
  Theta[n_groups]=1-sum(Theta[1:(n_groups-1)]); // to make a simplex.
}

model {
  real contributions[n_groups];
  real y_mu;
  real alpha = 1;

  // priors
  sigma ~ normal(.5,4);
  mu ~ normal(0,20);
  v ~ beta(1,alpha);
  for(i in 1:N){
    if(xcensl[i]==0 && xcensu[i]==0){
      xtrue[i] ~ normal(xobs[i]-0.5,2);
    }else if(xcensl[i]==1 && xcensu[i]==0){
      xtrue[i] ~ normal(xobs[i]-1.5,3);
    }else if(xcensl[i]==0 && xcensu[i]==1){
      xtrue[i] ~ normal(xobs[i]+.5,3);
    }
  }


	// update xtrue
	for(i in 1:N){
	    target+=xlikecens(xobs[i], xtrue[i],xsig,xcensu[i],xcensl[i]);
	}

  // update ytrue and logistic coefficients
  for(i in 1:N){
    y_mu = getytrue(coef,xtrue[i],yobs[i]);
    target+=ylikecens(yobs[i], y_mu, ysig, ycensu[i], ycensl[i]);
  }

  // update mixture normals
  for(i in 1:N) {
    for(k in 1:n_groups) {
      contributions[k] = log(Theta[k]) + normal_lpdf(xtrue[i] | mu[k], sigma[k]);
    }
    target += log_sum_exp(contributions);
  }

}

generated quantities{

  real MIC_Dens[Ngrid];
  real gx[Ngrid];
  real mb;
	real fx[Ngrid];

  // compute MIC density
  for (i in 1:Ngrid) {
	  MIC_Dens[i] = 0;
	  for (j in 1:n_groups) {
		  MIC_Dens[i]=MIC_Dens[i]+(Theta[j]*1/sqrt(2*pi()*sigma[j]^2)*
		    exp(-(xgrid[i]-mu[j])^2/(2*sigma[j]^2)));
	  }
  }

  //compute MIC-DIA relationship
	mb = (2*coef[3]*coef[4])/(coef[3]+coef[4]);
	for (i in 1:Ngrid){
  	fx[i] = 1/(1+exp(-mb*(coef[2]-xgrid[i])));
  	gx[i] = coef[1]*(fx[i]*exp(coef[3]*(coef[2]-xgrid[i]))+(1-fx[i])*exp(coef[4]*
  		(coef[2]-xgrid[i])))/(1+fx[i]*exp(coef[3]*(coef[2]-xgrid[i]))+(1-fx[i])*
  		exp(coef[4]*(coef[2]-xgrid[i])));
	}

}

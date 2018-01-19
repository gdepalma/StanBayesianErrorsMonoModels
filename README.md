# StanBayesianErrorsMonoModels

This was attempt to use STAN to estimate the models in the BayesianMonoErros package and described here: http://onlinelibrary.wiley.com/doi/10.1002/sim.7533/abstract

Implementation was somewhat successful.  FOr some data sets it worked well.  However, for other data sets the STAN model would not converge or take an extremely long time.  I believe this to be the number of unknown parameters (>1000) and the trying to approximate the Dirichlet Mixture of Normals distribution in STAN.

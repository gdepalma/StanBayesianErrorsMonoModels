# StanBayesianErrorsMonoModels

This was attempt to use STAN to for the models in the BayesianMonoErros package and described here: http://onlinelibrary.wiley.com/doi/10.1002/sim.7533/abstract

Implementation was somewhat successful.  However for some data sets the STAN would not converge or take an extremely long time.  I belive this to be the numebr of unknow parameters (>1000) and the trying to approximeate the Dirichlet Mixture of Normals distribution in STAN.

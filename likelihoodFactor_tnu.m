function dTC_tnu = likelihoodFactor_tnu(a,x,y,b)
nu = 1+1e6.*normcdf(a);
tau = 0.99.*(2.*normcdf(b)-1);

[dTC_tnu,~] = copulaDensities(x,y,tau,nu);
end
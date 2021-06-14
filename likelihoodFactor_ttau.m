function dTC_ttau = likelihoodFactor_ttau(a,x,y,b)
nu = 1+1e6.*normcdf(b);
tau = 0.99.*(2.*normcdf(a)-1);

[dTC_ttau,~] = copulaDensities(x,y,tau,nu);
end
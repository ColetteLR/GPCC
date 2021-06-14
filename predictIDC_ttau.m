function ttau = predictIDC_ttau(ret,zNew,nsamples)
[mZnew,vZnew] = predictFITC(ret.gFITCinfo, ret.f1Hat.eta2, ret.f1Hat.eta1, zNew);

rvs = randn(nsamples,1)';
nrvs = sqrt(vZnew).*rvs + mZnew;
ttau = (0.99 .* (2 .* normcdf(nrvs) - 1));
end

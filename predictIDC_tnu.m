function tnu = predictIDC_tnu(ret,zNew,nsamples)
[mZnew,vZnew] = predictFITC(ret.gFITCinfo, ret.f1Hat.eta2, ret.f1Hat.eta1, zNew);

rvs = randn(nsamples,1)';
nrvs = sqrt(vZnew).*rvs + mZnew;
tnu = (1 + 1e6 .* normcdf(nrvs));
end

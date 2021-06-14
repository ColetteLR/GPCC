function evalLL = evaluateLLCTC(ret, Xtest, Ytest, Ztest, nsamples)
if ~exist('nsamples','var')
    nsamples = 100000;
end

% Generate predictions for each datapoint

tau = predictIDC_ttau(ret.retTau, Ztest, nsamples);
nu = predictIDC_tnu(ret.retNu, Ztest, nsamples);
	
Xtest = matrix(Xtest,size(tau,1),size(tau,2));
Ytest = matrix(Ytest,size(tau,1),size(tau,2));

% Evaluate the test log likelihood on each data point
[dTC, ~] = copulaDensities(Xtest, Ytest, tau, nu);
evalLL.ll = log(dTC * (repmat(1 ./ size(tau,2),size(tau,2), 1)));

% Obtain the median and the quantiles for tau and nu
evalLL.meanTau = mean(tau,2);
evalLL.medianTau = median(tau,2);
evalLL.quantileInfTau = quantile(tau, 0.1, 2);
evalLL.quantileSupTau = quantile(tau, 0.9, 2);

evalLL.meanNu = mean(nu,2);
evalLL.medianNu = median(nu,2);
evalLL.quantileInfNu = quantile(nu, 0.1, 2);
evalLL.quantileSupNu = quantile(nu, 0.9, 2);

end

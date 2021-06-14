function [m,v] = predictFITC(gFITCinfo,tauTilde,muTilde,Xtest)

% Predictions using the generatised FITC approximation

% tauTilde = n-dimensional vector with the inverse variances  
%        of the Gaussian approximation to the likelihood factor
% muTilde = n-dimensional vector with the inverse variances times the mean
%        of the Gaussian approximation to the likelihood factor

% [m,v] = mean and variance of the predictions at the test points

[mm,dd] = size(Xtest);

Xtest = Xtest .* matrix(gFITCinfo.l,mm,dd,true); % by row
pStar = computeKnm(Xtest, gFITCinfo.XbarScaled, gFITCinfo.sigma);

dStar = computeDiagKnn(Xtest, gFITCinfo.sigma, gFITCinfo.sigma0) - ...
		((pStar*(gFITCinfo.R)').^2) * ones(gFITCinfo.m,1);

%Dnew = gFITCinfo.D ./ (1 + gFITCinfo.D .* tauTilde);
Pnew = matrix(1 ./ (1 + gFITCinfo.D .* tauTilde),gFITCinfo.n, gFITCinfo.m) .* ...
        gFITCinfo.P;

Rnew = (rot180( (chol(rot180(eye(gFITCinfo.m) + ...
		(gFITCinfo.PRt)'*(matrix(tauTilde ./ (1 + gFITCinfo.D .* tauTilde), ...
        gFITCinfo.n, gFITCinfo.m) .* gFITCinfo.PRt))))')) ...
        \ gFITCinfo.R; %backsolve for R

gammaNew = Rnew' * (Rnew * (Pnew' * (muTilde + gFITCinfo.AInvM0)));

% Obtain the new marginals

mPrediction = pStar * gammaNew;

% Add the contribution of the prior mean

mPrediction = mPrediction + gFITCinfo.m0 - pStar * gFITCinfo.RtRPAInvM0; % mean
vPrediction = dStar + (pStar * Rnew').^2 * ones(gFITCinfo.m,1); % variance

m = mPrediction;
v = vPrediction;
end

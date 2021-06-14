function NuTauLogZ = fitCTC(X, Z, m)
[T,~] = size(X);

if ~exist('m','var')
    m = round(T./10);
end

% Fit an unconditional copula to the data
[x,nlogZ] = fitUTC(X);
retOptimizerTau = x(1);
retOptimizerNu = x(2);

% Initalize the estimates of tau and nu

tau.m = repmat(norminv(((retOptimizerTau./0.99)+1)./2),T,1); 
tau.v = ones(T,1);

nu.m = repmat(log(retOptimizerNu-1),T,1); 
nu.v = ones(T,1);

% Store the logZ sequence
logZ = zeros(9,1); % preallocate for processing speed
logZ(1) = -nlogZ;

% Start the iterative algorithm

retTau = [];
retNu = [];
evidence = NaN;
for i=1:4

    % Estimate tau given the current estimate of nu

	retNu = epIDCexternal_tnu(X, Z, retOptimizerNu, tau, m, evidence, retNu);
    logZ(2*i) = retNu.logZ;
    
    nu.m = retNu.f2Hat.eta1 ./ retNu.f2Hat.eta2;
    nu.v = 1 ./ retNu.f2Hat.eta2;

	% Estimate nu given the current estimate of tau

	retTau = epIDCexternal_ttau(X, Z, retOptimizerTau, nu, m, retTau);
    logZ(2*i+1) = retTau.logZ;
    
    tau.m = retTau.f2Hat.eta1 ./ retTau.f2Hat.eta2;
    tau.v = 1 ./ retTau.f2Hat.eta2;
    
end

NuTauLogZ.retNu = retNu;
NuTauLogZ.retTau = retTau;
NuTauLogZ.logZsequence = logZ;

end

function [x,fval] = fitUTC(X)

theta0 = [0;10]; % initial parameter values
lb = [-0.99;1]; % parameter lower bound
ub = [0.99;1e6]; % parameter upper bound

% Optimise/ minimise the negative log-likelihood function
[x,fval] = fmincon(@neglogL,theta0,[],[],[],[],lb,ub); 

    % Calculate the negative log-likelihood
    function ll = neglogL(parameters)
    T=size(X,1);
    [~, logdTC] = copulaDensities(X(:,1),X(:,2),repmat(parameters(1),T,1),repmat(parameters(2),T,1));
    ll = -sum(logdTC);
    end

end

% Import dataset with uniformly distributed variables

X = importdata('udata.txt');
n = size(X,1);

% Adjust a conditional Student t copula to the first n-1 observations

Z = (1:n)';

Xtrain = X(1:size(X,1)-1,:);
Xtest = X(size(X,1),:);

Ztrain = Z(1:(size(Z,1)-1),:);
Ztest = Z(size(Z,1),:);

m=50;
ret_sim = fitCTC(Xtrain, Ztrain, m);

% Evaluate the test ll of the conditional copula

retTest_sim = evaluateLLCTC(ret_sim, Xtest(:,1), Xtest(:,2), Ztest);



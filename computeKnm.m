% Knm = nxm kernel matrix between the data points and pseudo-inputs
function Knm = computeKnm(X,Xbar,sigma)
n = size(X,1);
m = size(Xbar,1);
Q = matrix(sum(X.^2,2),n,m);
Qbar = matrix(sum(Xbar.^2,2),n,m,true); % by row
distance = Qbar+Q-2.*X*Xbar';
Knm = sigma.*exp(-0.5.*distance);
end
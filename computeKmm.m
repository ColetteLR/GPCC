% Kmm = mxm kernel (covariance) matrix for the pseudo-inputs
function Kmm = computeKmm(Xbar,sigma,sigma0)
m = size(Xbar,1);
Q = matrix(sum(Xbar.^2,2),m,m); % sum over rows
distance = Q+Q'-2.*Xbar*Xbar';
K = exp(-0.5.*distance);
Kmm = sigma.*K + sigma0.*eye(m) + 1e-4.*eye(m);
end
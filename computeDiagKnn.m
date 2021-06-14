% diagKnn = n-dimensional vector of the diagonal of the kernel matrix for data points
function diagKnn = computeDiagKnn(X,sigma,sigma0)
diagKnn = repmat(sigma+1e-4+sigma0,size(X,1),1);
end

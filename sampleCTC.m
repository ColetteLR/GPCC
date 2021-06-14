function sample = sampleCTC(n,ret,Ztest)
sample = zeros(n, 2);
tau = double(predictIDC_ttau(ret.retTau, matrix(Ztest, 1, length(Ztest)), n));
nu = double(predictIDC_tnu(ret.retNu, matrix(Ztest, 1, length(Ztest)), n));

for i=1:n
    sample(i,:) = copularnd('t', sin((tau(i) * pi) ./ 2), nu(i),1);
end
end
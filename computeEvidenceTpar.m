function logZret = computeEvidenceTpar(a,f1Hat,f2Hat,gFITCinfo)
logZ = sum(f1Hat.logZ,'all');

logZret = logZ + sum(0.5 .* log(2.*pi) - 0.5 .* log(f2Hat.eta2) + ...
    0.5 .* log(a.eta2) - 0.5 .*a.eta1.^2 ./ a.eta2 + ...
    0.5 .* f1Hat.eta1.^2 ./ f1Hat.eta2 + ...
    0.5 .* f2Hat.eta1.^2 ./ f2Hat.eta2,'all');

logZret = logZret + getFITCevidence(gFITCinfo,f1Hat.eta2,f1Hat.eta1);

end

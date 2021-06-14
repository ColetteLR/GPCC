% Density of Student's t copula in terms of Kendall's tau
function [dTC, logdTC] = copulaDensities(x,y,tau,nu)
rho=sin((tau.*pi)./2); 
x=tinv(x,nu); 
y=tinv(y,nu);

%dTC = copulapdf('t',u,rho,nu);
dTC = (1./(2.*pi.*sqrt(1-rho.^2))) .* ...
    (1+(x.^2+y.^2-2.*rho.*x.*y)./(nu.*(1-rho.^2))).^(-(nu+2)./2) .* ...
    (tpdf(x,nu).^(-1)).*(tpdf(y,nu).^(-1));

logdTC = -log(2*pi) - 0.5.*log(1-rho.^2) - ...
    (nu+2)./2 .* log(1+(x.^2+y.^2-2.*rho.*x.*y)./(nu.*(1-rho.^2))) - ...
    log(tpdf(x,nu)) - log(tpdf(y,nu));
end


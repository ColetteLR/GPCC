function PD = checkPDposterior(gFITCinfo,tauTilde)
M = rot180(eye(gFITCinfo.m)+gFITCinfo.PRt'* ...
    (matrix(tauTilde/(1+gFITCinfo.D.*tauTilde),gFITCinfo.n,gFITCinfo.m).*gFITCinfo.PRt));

eigenValues = eig(M);
eigenValues(abs(eigenValues)<1e-6) = 0;

if sum(any(eigenValues<=0),'all')>0
    PD=false;
    return
end

if sum(any((1+gFITCinfo.D .* tauTilde)<0),'all')>0
    PD=false;
    return
end

PD=true;
end

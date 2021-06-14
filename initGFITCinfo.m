% Initialise the function with the problem information
%function gFITCinfo = initGFITCinfo(X,Xbar,sigma,sigma0,l,m0)
function gFITCinfo = initGFITCinfo(X,Xbar,sigma,sigma0,l,m0)

if ~exist('m0','var')
    % parameter does not exist
    m0 = 0;
end

[m,d] = size(Xbar);
n = length(X);

% Initialise the structure with data and kernel hyper-parameter
gFITCinfo.X = X; % nxd data point
gFITCinfo.Xbar = Xbar; % mxd inducing (pseudo) inputs
gFITCinfo.m = m;
gFITCinfo.n = n;
gFITCinfo.d = d;
gFITCinfo.sigma = exp(sigma) + 1e-3; % log-amplitude of GP
gFITCinfo.sigma0 = exp(sigma0); % log-noise level in GP
gFITCinfo.l = sqrt(exp(l)); % d-dimensional log-lengthscales
gFITCinfo.m0 = m0; % mean of approximate prior

% Compute the kernel matrices
gFITCinfo.XbarScaled = gFITCinfo.Xbar .* matrix(gFITCinfo.l,m,d,true); % by row
gFITCinfo.XScaled = gFITCinfo.X .* matrix(gFITCinfo.l,n,d,true);
gFITCinfo.Kmm = computeKmm(gFITCinfo.XbarScaled,gFITCinfo.sigma,gFITCinfo.sigma0);
gFITCinfo.Knm = computeKnm(gFITCinfo.XScaled,gFITCinfo.XbarScaled,gFITCinfo.sigma);
gFITCinfo.P = gFITCinfo.Knm;    %P=Knm

M = gFITCinfo.Kmm;

%%% Error testing
if sum(isnan(M),'all')>0 || sum(isinf(M),'all')>0
   % error('matrix contains nan or inf values');
   disp('in')
   gFITCinfo = [];
   return
end

eigenValues = eig(M);
eigenValues(abs(eigenValues)<1e-6) = 0;

if sum(any(eigenValues<=0),'all')>0
    %error('negative eigen values');
    gFITCinfo = [];
    return;
end

%%%

gFITCinfo.R = cholInverse(gFITCinfo.Kmm);
gFITCinfo.PRt = gFITCinfo.P*gFITCinfo.R';

% Compute the diagonal matrices
gFITCinfo.diagKnn = computeDiagKnn(gFITCinfo.X,gFITCinfo.sigma,gFITCinfo.sigma0); 
gFITCinfo.D = gFITCinfo.diagKnn-((gFITCinfo.PRt).^2)*ones(m,1); % nx1 vector

% Compute (A^-1)m0 efficiently using the woodbury formula and cholesky
% representations
% M = Im + RP'(D^-1)PR'
M = eye(m)+gFITCinfo.PRt'*(matrix((gFITCinfo.D).^-1,n,m).*gFITCinfo.PRt);

%%% Error testing
if sum(isnan(M),'all')>0 || sum(isinf(M),'all')>0
    %error('matrix contains nan or inf values');
    gFITCinfo = [];
    return
end

eigenValues = eig(M);
eigenValues(abs(eigenValues)<1e-6) = 0;

if sum(any(eigenValues<=0),'all')>0
    %error('negative eigen values');
    gFITCinfo = [];
    return
end
%%%

% L0 = cholInv(Im + RP'(D^-1)PR')
L0 = cholInverse(M);

% A^-1 = (D^-1) - ((D^-1)PR')cholInv(Im + RP'(D^-1)PR')'cholInv(Im + RP'(D^-1)PR')(P'R(D^-1))
% (A^-1)m0 = (D^-1)m0 - ((D^-1)PR')cholInv(Im + RP'(D^-1)PR')'cholInv(Im + RP'(D^-1)PR')(P'R(D^-1))m0
gFITCinfo.AInvM0 = ((gFITCinfo.D.^-1).*gFITCinfo.m0 - ...
        (matrix(gFITCinfo.D.^-1,n,m).*gFITCinfo.PRt)* L0' * ...
         (L0*(gFITCinfo.PRt'*((gFITCinfo.D.^-1).*gFITCinfo.m0))));

gFITCinfo.RtRPAInvM0 = gFITCinfo.R'*(gFITCinfo.PRt'*gFITCinfo.AInvM0);

end

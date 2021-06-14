% Compute the derivative of the evidence with respect to the 
% kernel hyper-parameters and the pseudo inputs

% Outputs a list with the derivatives of the evidence
function dLogZd = ...
		 computeDerivativesEvidence(gFITCinfo, tauTilde, muTilde)

% Compute the auxiliary matrices and vectors

mHat = muTilde ./ tauTilde - gFITCinfo.m0;
A = (tauTilde.^-1 + gFITCinfo.D).^-1;
B = cholInverse(eye(gFITCinfo.m) + gFITCinfo.PRt' * ...
	(matrix(A,gFITCinfo.n, gFITCinfo.m) .* gFITCinfo.PRt));
C = matrix(A,gFITCinfo.n, gFITCinfo.m) .* (gFITCinfo.PRt * B');
e = (tauTilde.^-1 + gFITCinfo.D).^-1 .* mHat - C * (C' * mHat);
PRtR = gFITCinfo.PRt * gFITCinfo.R;
M1 = (PRtR' * C) * C' - PRtR' .* ...
	matrix(A,gFITCinfo.m, gFITCinfo.n,true);
M2 = 0.5 .* (M1 * PRtR);
v1 = (1 ./ sqrt(2)) .* (C.^2 * ones(gFITCinfo.m,1) - A);
v2 = (1 ./ sqrt(2)) .* ones(gFITCinfo.n,1);
v3 = (gFITCinfo.PRt * gFITCinfo.R)' * e;
e2 = (1 ./ sqrt(2)) .* e;
M1prima = M1 - 2 .* PRtR' .* ...
    matrix(v1 .* v2 + e2.^2, gFITCinfo.m, gFITCinfo.n,true); % by row
M2prima = (PRtR' .* ...
    matrix(v1 .* v2 + e2.^2, gFITCinfo.m, gFITCinfo.n,true)) * PRtR - M2;

% Evaluate the derivative of P, Kmm and diag(Knn) with respect to sigma

dPdSigma = gFITCinfo.P ./ gFITCinfo.sigma .* (gFITCinfo.sigma - 1e-3);
dKmmdSigma = (gFITCinfo.Kmm - gFITCinfo.sigma0 .* eye(gFITCinfo.m) - ...
	1e-4 .* eye(gFITCinfo.m)) ./ gFITCinfo.sigma .* ...
	(gFITCinfo.sigma - 1e-3);
dDiagKnndSigma = (gFITCinfo.diagKnn - gFITCinfo.sigma0 - 1e-4) ./ ...
	gFITCinfo.sigma .* (gFITCinfo.sigma - 1e-3);

% Compute the derivative of the log evience with respect to sigma

term1 = sum((v1 .* v2 + e2.^2 ) .* dDiagKnndSigma,"all");
term2 = sum(M1prima .* dPdSigma',"all");
term3 = sum(M2prima .* dKmmdSigma,"all");
term4 = double(e' * dPdSigma * v3);
term5 = double(-0.5 .* v3' * dKmmdSigma * v3);
dLogZdSigma = term1 + term2 + term3 + term4 + term5;

% Evaluate the derivative of P, Kmm and diag(Knn) with respect to sigma0

dPdSigma0 = zeros(gFITCinfo.n, gFITCinfo.m);
dKmmdSigma0 = gFITCinfo.sigma0 .* eye(gFITCinfo.m);
dDiagKnndSigma0 = gFITCinfo.sigma0 .* ones(gFITCinfo.n,1);

% Compute the derivative of the log evience with respect to sigma0

term1 = sum((v1 .* v2 + e2.^2 ) .* dDiagKnndSigma0,"all");
term2 = sum(M1prima .* dPdSigma0',"all");
term3 = sum(M2prima .* dKmmdSigma0,"all");
term4 = double(e' * dPdSigma0 * v3);
term5 = double(-0.5 .* v3' * dKmmdSigma0 * v3);
dLogZdSigma0 = term1 + term2 + term3 + term4 + term5;

% Compute the derivatives with respect to the lengthscales

dLogZdl = zeros(length(gFITCinfo.l),1);
for i=1:length(gFITCinfo.l)

	% Evaluate the derivative of P, Kmm and diag(Knn) with
	% respect to sigma0

	Q = matrix(gFITCinfo.XScaled(:,i).^2,gFITCinfo.n,gFITCinfo.m);
	Qbar = matrix(gFITCinfo.XbarScaled(:,i).^2,gFITCinfo.n,gFITCinfo.m,true);
	distance = Qbar + Q - 2 .* gFITCinfo.XScaled(:,i) * ...
		gFITCinfo.XbarScaled(:,i)';
	dPdl = -gFITCinfo.P .* distance;

    Q = matrix(gFITCinfo.XbarScaled(:,i).^2,gFITCinfo.m,gFITCinfo.m);
    Qbar = matrix(gFITCinfo.XbarScaled(:,i).^2,gFITCinfo.m,gFITCinfo.m,true);
    distance = Qbar + Q - 2 .* gFITCinfo.XbarScaled(:,i) * ...
    	gFITCinfo.XbarScaled(:,i)';
    dKmmdl = -(gFITCinfo.Kmm - ...
    	gFITCinfo.sigma0 .* eye(gFITCinfo.m) - ...
    	1e-4 .* eye(gFITCinfo.m)) .* distance;
    
    dDiagKnndl = zeros(gFITCinfo.n,1);

    % Compute the derivative of the log evience with respect to sigma0

	term1 = sum((v1 .* v2 + e2.^2 ) .* dDiagKnndl,"all");
	term2 = sum(M1prima .* dPdl',"all");
	term3 = sum(M2prima .* dKmmdl,"all");
	term4 = double(e' * dPdl * v3);
	term5 = double(-0.5 .* v3' * dKmmdl * v3);
	dLogZdl(i) = (term1 + term2 + term3 + term4 + term5) .* 0.5 .* ...
			gFITCinfo.l(i).^(-0.5);
end

% Compute the derivatives with respect to the pseudo inputs

dLogZdXbar = zeros(gFITCinfo.m, gFITCinfo.d);
for i=1:length(gFITCinfo.l)

	% Evaluate the derivative of P, Kmm and diag(Knn) with
	% respect to the pseudo-inputs

	distance = matrix(gFITCinfo.XScaled(:,i),gFITCinfo.n,gFITCinfo.m) - ...
        matrix(gFITCinfo.XbarScaled(:,i),gFITCinfo.n, gFITCinfo.m,true);
	dPdXbar = gFITCinfo.P .* distance .* gFITCinfo.l(i);

	distance = matrix(gFITCinfo.XbarScaled(:,i),gFITCinfo.m,gFITCinfo.m) - ...
        matrix(gFITCinfo.XbarScaled(:,i),gFITCinfo.m, gFITCinfo.m,true);
	dKmmdXbar = (gFITCinfo.Kmm - gFITCinfo.sigma0 .* eye(gFITCinfo.m) - ...
		1e-4 .* eye(gFITCinfo.m)) .* distance .* gFITCinfo.l(i);

	dDiagKnndXbar = zeros(gFITCinfo.n,1);

	% Compute the derivative of the log evience with respect to
	% the pseudo-inputs

	term1 = sum((v1 .* v2 + e2.^2) .* dDiagKnndXbar,"all") .* ones(gFITCinfo.m,1);
	term2 = (M1prima .* dPdXbar') * ones(gFITCinfo.n,1);
	term3 = -2 .* (M2prima .* dKmmdXbar) * ones(gFITCinfo.m,1);
	term4 = double(dPdXbar' * e .* v3);
	term5 = -dKmmdXbar' * v3 .* v3;
	dLogZdXbar(:,i) = term1 + term2 + term3 + term4 + term5;
end

% Compute the gradient with respect to the mean of the prior 
% Gaussian Process

dLogZdm0 = sum(e,"all");

dLogZd.dLogZdSigma = dLogZdSigma;
dLogZd.dLogZdSigma0 = dLogZdSigma0;
dLogZd.dLogZdl = dLogZdl;
dLogZd.dLogZdXbar = dLogZdXbar;
dLogZd.dLogZdm0 = dLogZdm0;

end

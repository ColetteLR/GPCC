function bestSoFar = epIDCexternal_ttau(X, Z, meanTau, nuValues, m, ret)

if ~exist('ret','var')
    % parameter does not exist
    ret = [];
end

m0 = norminv(((meanTau ./ 0.99) + 1) ./ 2);

% Initialize the hyper-parameters if no initial solution is given

if isempty(ret)
    sigma = 0;
    sigma0 = -10;
    estTtau = sigest(Z);
    l = log(sqrt(estTtau));
    Zbar = matrix(Z(randsample(size(Z,1),m),:),m,size(Z,2));
    eps = 0.05;
else
    sigma = ret.sigma;
    sigma0 = ret.sigma0;
    l = ret.l;
    Zbar = ret.Zbar;
    eps = ret.eps;
end

% Initialise the gradient optimisation process

try
    retNew = epIDCinternal_ttau(X, Z, Zbar, sigma, sigma0, l, m0, nuValues, ret);
catch
    bestSoFar = ret;
    return
end

ret = retNew;

bestSoFar = ret;
bestSoFar.eps = eps;

disp([num2str(0),' New evidence: ',num2str(ret.logZ)])

convergence = false;
iteration = 1;
while convergence == false && iteration < 50
    
    sigmaNew = sigma + eps .* ret.gradientLogZ.dLogZdSigma;
	sigma0New = sigma0 + eps .* ret.gradientLogZ.dLogZdSigma0;
	lNew = l + eps .* ret.gradientLogZ.dLogZdl;
	ZbarNew = Zbar + eps .* ret.gradientLogZ.dLogZdXbar;

    try
        gFITCinfo = initGFITCinfo(Z, ZbarNew, sigmaNew, sigma0New, lNew, m0);
    catch
        %bestSoFar = bestSoFar;
        return
    end
	
    counter = 1;
	while isempty(gFITCinfo) || checkPDposterior(gFITCinfo, ret.f1Hat.eta2) == false
		eps = eps .* 0.5;
		sigmaNew = sigma + eps .* ret.gradientLogZ.dLogZdSigma;
		sigma0New = sigma0 + eps .* ret.gradientLogZ.dLogZdSigma0;
		lNew = l + eps .* ret.gradientLogZ.dLogZdl;
		ZbarNew = Zbar + eps .* ret.gradientLogZ.dLogZdXbar;

		try
            gFITCinfo = initGFITCinfo(Z, ZbarNew, sigmaNew, sigma0New, lNew, m0);
        catch
            %bestSoFar = bestSoFar;
            return
		end

		if counter > 1e2
            %bestSoFar = bestSoFar;
			return
		end
        
        counter = counter + 1;
	end
    
    sigma = sigmaNew;
    sigma0 = sigma0New;
	l = lNew;
	Zbar = ZbarNew;

    try
        retNew = epIDCinternal_ttau(X, Z, Zbar, sigma, sigma0, l, m0, nuValues, ret);
    catch
        %bestSoFar = bestSoFar;
        return
    end

	if (abs(retNew.logZ - ret.logZ) / abs(ret.logZ)) < 1e-4
		convergence = true;
        %break
	end
    
    if (retNew.logZ < ret.logZ)
    	eps = eps .* 0.5;
    else
    	eps = eps .* 1.1;
    end
    
    if (retNew.logZ > bestSoFar.logZ)
    	bestSoFar = retNew;
    	bestSoFar.eps = eps;
    end
    
    disp([num2str(iteration), ' New evidence: ', num2str(retNew.logZ), ' eps: ', num2str(eps)])
    	
    ret = retNew;
    
    iteration = iteration + 1;
    
end

%bestSoFar=bestSoFar;
end


 %function [a, f1Hat, f2Hat, logZ, gradientLogZ, gFITCinfo, ...
%		sigma, sigma0, l, m0, Zbar] = ...
function internal_ttau = epIDCinternal_ttau(X, Z, Zbar, sigma, sigma0, l, m0, nuValues, start)

if ~exist('start','var')
    % parameter does not exist
    start = [];
end

% Initialize the structure with the problem information

gFITCinfo = initGFITCinfo(Z, Zbar, sigma, sigma0, l, m0);

% Initialize the approximate factors

f1Hat.eta1 = repmat(1e-100,gFITCinfo.n,1); 
f1Hat.eta2 = repmat(1e-100,gFITCinfo.n,1); 
f1Hat.logZ = zeros(gFITCinfo.n,1);

% Initialize the posterior approximation

%a.eta1 = zeros(gFITCinfo.n,1); a.eta2 = zeros(gFITCinfo.n,1);

% Refine the second approximate factor

f2Hat.eta2 = gFITCinfo.diagKnn.^(-1);
a.eta2 = f2Hat.eta2;
f2Hat.eta1 = m0 ./ gFITCinfo.diagKnn;
a.eta1 = f2Hat.eta1;

% Check for an initial solution

if ~isempty(start)
    a = start.a;
    f1Hat = start.f1Hat;
    f2Hat = start.f2Hat;
end

% Main loop of EP

i = 1;
damping = 1;
convergence = false;
while convergence == false && i<50
    aOld = a;
    
    % Refine the first approximate factor
    
    f1HatOld = f1Hat;
    
    PD = false; % not positive definite
    while PD == false
        f1Hat = f1HatOld;
        
        % Create the auxiliary list m
        
        m.X = X(:,1);    
        m.Y = X(:,2);
        m.ma = f2Hat.eta1 ./ f2Hat.eta2;     
        m.va = f2Hat.eta2.^(-1);
        m.mb = nuValues.m;
        m.vb = nuValues.v;
		m.f1HatEta1 = f1Hat.eta1;    
        m.f1HatEta2 = f1Hat.eta2;     
        m.n = gFITCinfo.n;  
        m.logZ = f1Hat.logZ;
        m.damping = damping;
        
        m = refineLikelihoodFastTau01(m);
        
        f1Hat.eta1 = m.f1HatEta1;
        f1Hat.eta2 = m.f1HatEta2;
        f1Hat.logZ = m.logZ;

        PD = checkPDposterior(gFITCinfo,f1Hat.eta2);
        if PD == true
            break
        else
            damping = damping .* 0.5;
        end
        
    end  
    
    % Refine the second approximate factor
    
    [mNew, vNew] = computeTitledDistribution(gFITCinfo,f1Hat.eta2,f1Hat.eta1);
    eta1HatNew = mNew ./ vNew - f1Hat.eta1;
    eta2HatNew = vNew.^(-1) - f1Hat.eta2;
    
    index = find(~isfinite(eta2HatNew));
    if ~isempty(index)
        eta2HatNew(index) = f2Hat.eta2(index);
        eta1HatNew(index) = f2Hat.eta1(index);
    end    
    
    index = find(~isfinite(eta1HatNew));
    if ~isempty(index)
        eta2HatNew(index) = f2Hat.eta2(index);
        eta1HatNew(index) = f2Hat.eta1(index);
    end  
    
    index = find((damping .* eta2HatNew + (1 - damping) .* f2Hat.eta2 + f1Hat.eta2) < 0);
    if ~isempty(index)
        eta2HatNew(index) = f2Hat.eta2(index);
        eta1HatNew(index) = f2Hat.eta1(index);
    end 
    
    eta1HatNew(eta2HatNew <= 0) = f2Hat.eta1(eta2HatNew <= 0);
    eta2HatNew(eta2HatNew <= 0) = f2Hat.eta2(eta2HatNew <= 0);
    
    f2Hat.eta1 = damping .* eta1HatNew + (1 - damping) .* f2Hat.eta1;
    f2Hat.eta2 = damping .* eta2HatNew + (1 - damping) .* f2Hat.eta2;
    
    % Update the posterior approximation
    
    a.eta1 = f1Hat.eta1 + f2Hat.eta1;
    a.eta2 = f1Hat.eta2 + f2Hat.eta2;
    
    % Check for convergence
    
    change = max(abs(aOld.eta1 ./ aOld.eta2 - a.eta1 ./ a.eta2));
	change = max(change, max(abs(aOld.eta2.^(-1) - a.eta2.^(-1))));
	if change < 1e-4
        convergence = true;
        %break
    end

	% Annealed damping scheme
	damping = damping .* 0.95;

	i = i + 1;
end

% Compute the evidence and its gradient

logZ = computeEvidenceTpar(a, f1Hat, f2Hat, gFITCinfo);
gradientLogZ = computeDerivativesEvidence(gFITCinfo, f1Hat.eta2, f1Hat.eta1);

internal_ttau.a=a; 
internal_ttau.f1Hat=f1Hat; 
internal_ttau.f2Hat=f2Hat;
internal_ttau.logZ=logZ;
internal_ttau.gradientLogZ=gradientLogZ;
internal_ttau.gFITCinfo=gFITCinfo;
internal_ttau.sigma=sigma;
internal_ttau.sigma0=sigma0;
internal_ttau.l=l;
internal_ttau.m0=m0;
internal_ttau.Zbar=Zbar;

end


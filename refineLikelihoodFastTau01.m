function LTau = refineLikelihoodFastTau01(m)
% Auxiliary variables:

%int i;
%double Z, m1, m2;
%double maNew, vaNew, eta1New, eta2New, eta1HatNew, eta2HatNew;

X=m.X;    Y=m.Y;
ma=m.ma;    va=m.va;
mb=m.mb;    vb=m.vb;
f1HatEta1=m.f1HatEta1;    
f1HatEta2=m.f1HatEta2;
damping=m.damping;     
n=m.n;  
logZ=m.logZ;

%natural parameters
for i=1:n

    if (va(i) > 0 && vb(i) > 0)

		% Compute the normalization contant
		[Z, m1, m2] = computeNormalizationConstantAndMomentsTau(X(i), Y(i), mb(i), ma(i), va(i), 10);

		if (Z > 1e-10) 
            
			maNew = m1 ./ Z;
			vaNew = m2 ./ Z - maNew .* maNew;

			eta1New = maNew ./ vaNew; %shift = mean/variance
			eta2New = 1 ./ vaNew; %precision = inv(variance)
			eta1HatNew = eta1New - ma(i) ./ va(i); % cavity
			eta2HatNew = eta2New - 1 ./ va(i);

			if (eta2HatNew < 1e10 && eta1HatNew < 1e10 && eta1HatNew > -1e10)
				f1HatEta2(i) = damping .* eta2HatNew + (1 - damping) .* f1HatEta2(i); %update
				f1HatEta1(i) = damping .* eta1HatNew + (1 - damping) .* f1HatEta1(i);
				logZ(i) = log(Z);
			end
		end	
    end
end

LTau.f1HatEta2 = f1HatEta2;
LTau.f1HatEta1 = f1HatEta1;
LTau.logZ = logZ;
end

% Compute the normalization constant (moments)

function [Z, m1, m2] = computeNormalizationConstantAndMomentsTau(X, Y, a, ...
    mb, vb, nSplits)

    sqrtvb = sqrt(vb);

	maxb = mb + 6 .* sqrtvb;
	minb = mb - 6 .* sqrtvb;
	hb = (maxb - minb) ./ (nSplits .* 4);
    
    % Boole's rule for integration
    % (2*h/45)*(7 32 12 32 7)*(function values)
    
    boole = [7; repmat([32 12 32 7.*2]',nSplits-1,1); [32 12 32 7]'];
    weight = boole .* hb .* 2.0 ./ 45;
    
    xb = [minb; repmat(hb,length(weight)-1,1)];
    xb = cumsum(xb);
        
	densityValue = tDensityTau(X,Y,xb,repmat(a,length(xb),1)) .* ...
                normpdf(xb, mb, sqrtvb);
        
    Z = weight .* densityValue;
    m1 = weight .* densityValue .* xb;
    m2 = weight .* densityValue .* xb .* xb;
        
    Z=sum(Z,'all');
    m1=sum(m1,'all');
    m2=sum(m2,'all');
    
end

function dTC_ttau = tDensityTau(x,y,a,b)
nu = 1+1e6.*normcdf(b); % area under normal curve
tau = 0.99.*(2.*normcdf(a)-1);

[dTC_ttau,~] = copulaDensities(x,y,tau,nu);
end
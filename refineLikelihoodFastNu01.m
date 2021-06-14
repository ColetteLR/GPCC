function LNu = refineLikelihoodFastNu01(m)

% Auxiliary variables:

X=m.X;    Y=m.Y;
ma=m.ma;    va=m.va;
mb=m.mb;    vb=m.vb;
f1HatEta1=m.f1HatEta1;    
f1HatEta2=m.f1HatEta2;
damping=m.damping;     
n=m.n;  
logZ=m.logZ;

for i=1:n
    
    if (va(i) > 0 && vb(i) > 0)

		% Compute the normalization contant
        %[Z, m1, m2] = computeNormalizationConstantAndMomentsNuExternal(X(i), Y(i), ma(i), va(i), mb(i), vb(i), 8);
		[Z, m1, m2] = computeNormalizationConstantAndMomentsNuInternal(X(i), Y(i), ma(i), mb(i), vb(i), 10);
       
		if (Z > 1e-10) 

			maNew = m1 ./ Z;
			vaNew = m2 ./ Z - maNew .* maNew;

			eta1New = maNew ./ vaNew;
			eta2New = 1 ./ vaNew;
			eta1HatNew = eta1New - mb(i) ./ vb(i);
			eta2HatNew = eta2New - 1 ./ vb(i);

			if (eta2HatNew < 1e10 && eta1HatNew < 1e10 && eta1HatNew > -1e10)
				f1HatEta2(i) = damping .* eta2HatNew + (1 - damping) .* f1HatEta2(i);
				f1HatEta1(i) = damping .* eta1HatNew + (1 - damping) .* f1HatEta1(i);
				logZ(i) = log(Z);
			end
		end	
    end
end

LNu.f1HatEta2 = f1HatEta2;
LNu.f1HatEta1 = f1HatEta1;
LNu.logZ = logZ;

end


% Compute the normalization constant

function [Z, m1, m2] = computeNormalizationConstantAndMomentsNuInternal(X, Y, a, ...
    mb, vb, nSplits)

    sqrtvb = sqrt(vb);

	maxb = mb + 6 .* sqrtvb;
	minb = mb - 6 .* sqrtvb;
	hb = (maxb - minb) ./ (nSplits .* 4);
    
    % Boole's rule for integration:
    % (2*h/45)*(7 32 12 32 7)*(function values)
    
    boole = [7; repmat([32 12 32 7.*2]',nSplits-1,1); [32 12 32 7]'];
    weight = boole .* hb .* 2.0 ./ 45;
    
    xb = [minb; repmat(hb,length(weight)-1,1)];
    xb = cumsum(xb);
        
	densityValue = tDensityNu(X,Y,repmat(a,length(xb),1),xb) .* ...
                normpdf(xb, mb, sqrtvb);
            
    Z = weight .* densityValue;
    m1 = weight .* densityValue .* xb;
    m2 = weight .* densityValue .* xb .* xb;
        
    Z=sum(Z,'all');
    m1=sum(m1,'all');
    m2=sum(m2,'all');
    
end

function dTC_tnu = tDensityNu(x,y,a,b)
nu = 1 + 1e6 .* normcdf(b); %a
tau = 0.99 .* (2 .* normcdf(a) - 1); %b

[dTC_tnu,~] = copulaDensities(x,y,tau,nu);
end

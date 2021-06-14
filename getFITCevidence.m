% Evaluate the normalization constant of the product
% of the FITC prior and a multivariate Gaussian density with diagonal
% correlation matrix

% logZ = Logarithm of the normalization constant
function logZ = getFITCevidence(gFITCinfo, tauTilde, muTilde)

Dnew = gFITCinfo.D ./ (1 + gFITCinfo.D .* tauTilde);
Pnew = matrix(1 ./ (1 + gFITCinfo.D .* tauTilde),gFITCinfo.n,gFITCinfo.m) ...
        .* gFITCinfo.P;

% Rnew = (rot180((chol(rot180(eye(gFITCinfo.m) + (gFITCinfo.PRt)' * ...
%         (matrix(tauTilde ./ (1 + gFITCinfo.D .* tauTilde),gFITCinfo.n, gFITCinfo.m) .* gFITCinfo.PRt))))')) ...
% 		\gFITCinfo.R;
    
S = (rot180(eye(gFITCinfo.m) + (gFITCinfo.PRt)' * ...
    (matrix(tauTilde ./ (1 + gFITCinfo.D .* tauTilde),gFITCinfo.n, gFITCinfo.m) .* ...
     gFITCinfo.PRt)));

[R,flag] = chol(S);

if ~flag
    Rnew = (rot180(R'))\gFITCinfo.R;
    
    %disp('Factorizations successful.')
else
    [V,D] = eig(S);       % Calculate the eigendecomposition of your matrix (A = V*D*V')
    % where "D" is a diagonal matrix holding the eigenvalues of your matrix "A"
    d= diag(D);           % Get the eigenvalues in a vector "d"
    d(d <= 1e-7) = 1e-7;  % Set any eigenvalues that are lower than threshold "TH" ("TH" here being
    % equal to 1e-7) to a fixed non-zero "small" value (here assumed equal to 1e-7)
    D_c = diag(d);        % Built the "corrected" diagonal matrix "D_c"
    S_PD = V*D_c*V';      % Recalculate your matrix "A" in its PD variant "A_PD"
    Rnew = (rot180(S_PD'))\gFITCinfo.R;
    
    %disp('Factorizations failed.')
end

upsilon = gFITCinfo.AInvM0 + muTilde;
aNew = Dnew .* upsilon;
gammaNew = Rnew' * (Rnew * (Pnew' * upsilon));

% Obtain the new marginal means

mNew = double(aNew) + double(Pnew * gammaNew);

% (Rnew appears as a consequence of the matrix determinant lemma)

logZ = -gFITCinfo.n ./ 2 .* log(2 .* pi) + sum(log(diag(Rnew)),'all') - ...
		sum(log(diag(gFITCinfo.R)),'all') - 0.5 .* sum(log(1 + tauTilde .* ...
		gFITCinfo.D),'all') + 0.5 .* sum((muTilde + gFITCinfo.AInvM0) .* mNew,'all') - ...
        0.5 .* sum(muTilde.^2 ./ tauTilde,'all') - 0.5 .* ...
        sum(gFITCinfo.m0 .* gFITCinfo.AInvM0,'all');

end

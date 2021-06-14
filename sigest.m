function tPar = sigest(x, frac)
if ~exist('frac','var')
    % parameter does not exist
    frac = 0.5;
end

m = size(x,1);
n = floor(frac.*m);
index = randsample(m,n,true); % sampling with replacement
index2 = randsample(m,n,true);
temp = (x(index,:) - x(index2,:));
dist = sum((temp.^2),2);
srange = 1./quantile(dist(dist~=0),[0.9,0.5,0.1]);
tPar = mean(srange([1,3]));
end
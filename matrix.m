function mat = matrix(X,n,m,byrow)
if ~exist('byrow','var') 
    byrow = false; % set default
end

% Reproduce the R matrix function to ensure identical output

% The matrix can be shaped from X by row (byrow=true) or 
% column (byrow=false).

% Since the number of elements must not change when using reshape,
% create a vector of appropriate length using repmat.

% This vector does not have to be a sub-multiple or multiple of the number
% of rows in X, but a warning message is displayed if this is the case.


[T,d] = size(X);
if byrow == false
    % stack columns to form row vector
    X1 = reshape(X,T*d,1); % stack columns to form row vector
    % create a vector of appropriate length
    vector = repmat(X1,ceil((n*m)./length(X1)),1);
    mat = reshape(vector(1:n*m,:), n, m);
    
end
if byrow == true
    X1 = reshape(X,T*d,1);
    vector = repmat(X1,ceil((n*m)./length(X1)),1);
    mat = reshape(vector(1:n*m,:), m, n)';
end

% Warning message

%if mod((T*d),n) ~=0
%    warning(['In matrix(X, ',num2str(n),', ', num2str(m),'): ',...
%       'data length (',num2str(T*d),') is not a sub-multiple or multiple of the number of rows (', num2str(n),').'])
%end
end
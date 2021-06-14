% Cholesky decomposition of the inverse of a particular matrix
% @param	A	m x m positive definite matrix
% @return	L	m x m upper triangular matrix such that M^-1 = L*L'

function cholInv = cholInverse(A)
cholInv = rot180(chol(rot180(A))')\eye(length(A)); % chol(inv(A))
end

function Q = constructCoeffMatrixC(P1, P2, W)
% P1: points in view 1, each column is a point
% P2: points in view 2, each column is a point
% W: weights for each correspondence
%    --- 

% #observations
n_obs = size(P2, 2);

if nargin < 3
    % set all the weights to one
    W = ones(n_obs, 1);
end

% check dimension compatibility
if (n_obs < 1) || (n_obs ~= size(P1, 2))
    fprintf('ERROR: dimension of the set of points are NOT equal\n');
    return
end


Q = zeros(9);

for ii=1:n_obs
   % retrieve observations
   f = P1(:, ii);
   f_p=P2(:, ii);
   % compute kronecker and weight
   kronf = kron(f_p, f);
   C_i = symmetrize(W(ii) * (kronf) * (kronf)');
   % Q_{e, e}. 
   Q(1:9, 1:9) =  Q(1:9, 1:9) + C_i;
end

Q = symmetrize(Q);
Q = Q ./ trace(Q);
end %end-function

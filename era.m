function [Ard,Brd,Crd] = era(U,S,V,r,m,p)
%% This is the Eigensystem Realization Algorithm in 
% S.-Y.Kung, 1979, "A new identification and model reduction algorithm via singular
% value decomposition" 
% The difference is, that this is a general formulation, and does not
% require U*S*V = svd(Hankel)
%
%  INPUTS
%
%       U   matrix with orthogonal columns
%       S   diagonal matrix with "scaling" factors (i.e. sing. values)
%       V   matrix with orthogonal columns
%       r   desired reduced model order
%       m   input dimension of the system
%       p   output dimension of the system
%
%  OUTPUTS
%
%     Ard   Disrete time system matrix
%     Brd   disrete time input matrix
%     Crd   discrete time output matrix
%
% Copyright (c) MIT, 2016
% Boris Kramer(bokr@vt.edu), first draft written at Virginia Tech, 2013
%% ------------------------------------------------------------------------

nrH = max(size(U));
ncH = max(size(V));

S = diag(S);

% Truncate and partition orthogonal matrices U,V
U1 = U(1:nrH-p, 1:r);   % Cutting last block-row
U2 = U(p+1:nrH, 1:r);   % Cutting first block-row
V1 = V(1:ncH-m , 1:r);

sqS    = sqrt(S(1:r));
invsqS = ones(r,1)./sqS;

Ubar = U1'*U2;
Ard = Ubar.*(invsqS*sqS');
Brd = (sqS*ones(1,m)).*V1(1:m,:)';  
Crd = U1(1:p,:).*(ones(p,1)*sqS');  

end

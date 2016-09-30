function [I,J,U,S,V,Hdiff_full,condH,Hdiff_iterates,HnormIterates] = crossApprox(H,r,relTolH,epsTol,compMore)
%% Computes the cross approximation of a matrix via maximum volume concept
% The result will be:  indices I,j, such that:  H \approx H(I,:)*H(I,J)^(-1)*H(:,J) 
%
%    INPUTS
%
%       H               Hankel matrix (or more generally any LR matrix
%       r               rank estimate of H
%       relTolH         relative tolerance for ||H_{k} - H_{k-1}||\||H_{k}|| at iteration k		
%       compMore        computes more metrics, such ||H - H_{k}||, vol(H_{k}), cond(H_{k})
%       epsTol          tolerance for maxvol algorithm, typically 2e-2
%
%   OUTPUTS
%	
%       I,J             row and column index of H(I,J)
%       U,S,V           svd(H(I,J)) = U*S*V', needed for extra computations, and to define Gramians
%       Hdiff_iterates	||H_{k} - H_{k-1}||\||H_{k}|| for each iteration of the algorithm
%       HnormIterates	||H_{k}||_F for each iteration
%       condH           vector of condition numbers of H_{k} at every it.
%       Hdiff_full      ||H - H_{k}||\||H_{k}||
%
%
%  Copyright (c) MIT, 2016
%  Boris Kramer (bokramer@mit.edu) and Alex A. Gorodesky(goroda@mit.edu)
%% ------------------------------------------------------------------------

verbose = 0; % Set to 1 if comments needed

if compMore == 1  % get some additional quantitative information (expensive)
    Honorm = norm(H,'fro');
end

% Initializing output quantities
[nr,nc] = size(H);
J = randperm(nc,r); 
converged = 0;
iteration = 0;
Hdiff_full = [];
Hdiff_iterates =[];
Hdiff_iterates(1) = 0;
HnormIterates =[];
vol =[];
condH =[];


while converged == 0
    
    iteration = iteration + 1;
    
    [Q1,R1] = qr(H(:,J),0);
    I = maxvol(Q1, epsTol, verbose);
    [Q2,R2] = qr(H(I,:)',0);
    J = maxvol(Q2,epsTol, verbose);
    %I = J; % for symmetry preserving methods
    
    [U,S,V] = svd(H(I,J));
    Hinv =  V*(diag(1./diag(S)))*U';
    
    Cnew = H(:,J)*Hinv;
    Rnew = H(I,:);
    
    G1 = zeros(1,r^2); G2 = zeros(r^2,1);
    for ii = 1:nr
        G1 = G1 + kron(Cnew(ii,:),Cnew(ii,:));
        G2 = G2 + kron(Rnew(:,ii),Rnew(:,ii));
    end
    
    HnormIterates(iteration) = sqrt(G1*G2);   % same as norm(Hnew,'fro'), but cheaper
    
    % Fast computation of the difference of two LR matrices, see appendix A
    if iteration > 1 
        G1 = zeros(1,4*r^2); G2 = zeros(4*r^2,1);
        C = [Cnew, -Cold];
        R = [Rnew; Rold];
        for ii = 1:nr
            G1 = G1 + kron(C(ii,:),C(ii,:));
            G2 = G2 + kron(R(:,ii),R(:,ii));
        end

        Hdiff_iterates(iteration) = real(sqrt(G1*G2))/HnormIterates(iteration);
    end
    
    if (Hdiff_iterates(iteration) <  relTolH)  && (iteration > 1)
            converged = 1;
    end
    
   Cold = Cnew; 
   Rold = Rnew;
   
   if compMore == 1 % If we want more outputs 
     	Hnew = H(:,J)*Hinv*H(I,:);
        Hdiff_full(iteration) = norm(Hnew-H,'fro')/Honorm;
        vol(iteration) = log(abs(det(H(I,J))));    
        condH(iteration) = cond(R1*Hinv*R2');
        Hold = Hnew;
        Honorm = norm(Hold,'fro');
   end
   
   display(strcat('iteration ', num2str(iteration), ':  Difference__', num2str(Hdiff_iterates(iteration)) ) ); 
   
end

end


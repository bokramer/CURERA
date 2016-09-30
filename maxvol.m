function I = maxvol(Q,epsTol,verbose)
%% Computes the maximum volume approximation of the matrix Q
%
%	INPUTS 	
%
%       Q    	Matrix to be approximated
%       epsTol	tolerance to find pivot element (typically 2e-2)
%       verbose	set to 1 if you want infos printed out every iteration
%		
%   OUTPUTS 	
%	
%           I	column index for the maxvol approximation 
%
%
%   NOTES
%      
%   Based on the paper TT-cross approximation for multidimensional arrays (2010)
%   by Ivan Oseledets and Eugene Tyrtyshnikov 
%
%  Copyright (c) MIT, 2016
%  Boris Kramer (bokramer@mit.edu) and Alex A. Gorodesky(goroda@mit.edu)
%% -------------------------------------------------------------

if (nargin == 2)
    verbose = 0;
end

Qtt = Q;
[n,r] = size(Q);
eyemat = speye(r);
eyemat2 = speye(n);

[~,~,p] = lu(Q,'vector');

Qt = Qtt(p,:);
Qs = Qt(1:r,:);
[u,s,v] = svd(Qs);
Qsinv =  v*(diag(1./diag(s)))*u';
B = Qt*Qsinv;

converged = 0;
iteration = 0;
while converged == 0
    iteration = iteration + 1;
    
    absB = abs(B);
    [maxEachCol,maxElrow] = max(absB);
    [maxElvalue, maxElcol] = max(maxEachCol);
    maxElrow = maxElrow(maxElcol);
    if (verbose == 1)
        display(strcat('MaxVol Iteration:', num2str(iteration), 'Maximum element in B is: ', num2str(maxElvalue)))
    end
    if maxElvalue > (1+epsTol)
        tempCol = p(maxElrow);
        p(maxElrow) = p(maxElcol);
        p(maxElcol) = tempCol;
        
        Btemp = (B(r+1:end,maxElcol)+eyemat2(r+1:end,maxElrow))*(B(maxElrow,:)-eyemat(maxElcol,:));
        B(r+1:end,:) = B(r+1:end,:) - Btemp/B(maxElrow,maxElcol);
        
    else
        converged = 1;
    end
    
    if iteration > 10*r   % Automatic return if too many iterations needed
        converged = 1;    % matrix probably not low rank
    end
end

I = p(1:r);

end

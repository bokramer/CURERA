function Z = block_Hankel(C, R)
%% This function assembles the Hankel matrix from a series of Markov parameters
%
%   INPUTS
%
%        C      cell array of (n) ma \times na matrices - the first column of Z
%        R      cell array of (n) ma \times na matrices - the last row of Z
%
%
%   OUTPUTS 
%        Z      n \times n block Hankel, constant (blocks) along the anti-diagonals
%
% NOTE: we use C{1:n} and R{2:n}, the 'repeat' R{1} = C{n} is ignored
%
%
% Copyright (c) Virginia Tech, 2014
% Eugene Cliff and Boris Kramer (bokr@vt.edu), 2014
%% Gather size information and pre-allocate
   [mc, nc] = size(C{1});   
   ic = 1:mc; jc = 1:nc; 
   lenC = length(C);
   
   
   if nargin == 2 && lenC ~= length(R)
       error('The length(s) of arrays  C  and R must be equal \n');
   end
    
   Z = zeros(mc*lenC, nc*lenC); % pre-allocate
 
%% Fill in the upper anti-diagonal blocks 
%  upward from the first column to the first row

   for ii=1:lenC
       if any( size(C{ii}) ~= [ mc nc] )
           error('The C{%5i} array is not compatible with the C{1}\n', ii);
       end
       iz = ic + ii*mc;
       jz = jc - nc;
       for kk=1:ii
           iz = iz - mc;
           jz = jz + nc;
           Z(iz, jz) = C{ii};
       end
   end
   
   if nargin == 2
       
%% Fill in the lower anti-diagonal blocks
%  downward from last the column to the last row

     for ii=2:lenC
         if any( size(R{ii}) ~= [ mc nc] )
           error('The R{%5i} array is not compatible with the C{1}\n', ii);
         end
         iz = ic + (ii-2)*mc;
         jz = jc +  lenC*nc;
         for kk=ii:lenC
           iz = iz + mc;
           jz = jz - nc;
           Z(iz, jz) = R{ii};
         end
     end
   
   end

end


%% To see what this file does, do 

% C=cell(1,4);
%     for i=1:4 
%         C{i} = i;
%     end
%     
% C = 
% 
%     [1]    [2]    [3]    [4]
% R = 
%
%     [4]    [5]    [6]    [7] % Note, need the "Repeat" value R{1}=4;
% 
% Z = block_Hankel(C)
% 
% Z =
% 
%      1     2     3     4
%      2     3     4     0
%      3     4     0     0
%      4     0     0     0
% 
% Z = block_Hankel(C,R)
% 
% Z =
% 
%  1     2     3     4
%  2     3     4     5
%  3     4     5     6
%  4     5     6     7


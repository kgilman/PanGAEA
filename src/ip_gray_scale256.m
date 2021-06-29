function B=gray_scale256(A)
% B=gray_scale256(A)
% scales image array A to an 8-bit gray scale
% Does not truncate.
% It shifts by most negative value then compresses
% based upon the new max value to 256 levels.

% get dimensions of A
if ndims(A)>2
    error('Error: input cannot have more than 2 dimensions.')
end

L=256; % number of gray levels
L1=L-1;
A=double(A); %just in case

Amin=min(min(A));

B=A-Amin; % shift by the minimum value
Bmax=max(max(B));
scale=L1/Bmax;  % compression or expansion factor
B=B*scale; % force to L-1 max value
B=round(B);  % force to integer--is "floor" better? Can B ever exceed L1?
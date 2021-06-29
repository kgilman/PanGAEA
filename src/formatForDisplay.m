function [L, E, S] = formatForDisplay(L,E,S,M)
%
% Prepares a decomposition for display by doing the following:
%   - Clamping the background to [0, 1]
%   - Converting the corruption and sparse components to [0, 1], where 0.5
%     (gray) corresponds to no signal
%
% L: (ny x nx x nt) background
% E: (ny x nx x nt) corruption
% S: (ny x nx x nt) foreground
% M: (ny x nx x nt) frame masks
%

% Format background
L = max(0,min(L,1));
ML = repmat(any(M,3),[1, 1, size(M,3)]);
L(~ML) = 0;

% Format corruption
E = max(0,min(E + 0.5,1));
E(~M) = 0;

% Format foreground
S = max(0,min(S + 0.5,1));
% S = min(abs(S),0.5) * 2; %this is to be left uncommented
S(~M) = 0;

% L(M) = L(M) - min(L(M),[],1);
% L(M) = 255.*(L(M) ./ max(L(M),1));
% 
% E(M) = E(M) - min(E(M),[],1);
% E(M) = 255.*(E(M) ./ max(E(M),1));
% 
% S(M) = S(M) - min(S(M),[],1);
% S(M) = 255.*(S(M) ./ max(S(M),1));



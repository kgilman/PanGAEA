function [ e, s, xi, w, y1,y2,y3,h1,h2,h3, iter] = admm_srp_pangaea( U, v, Bf, C, idx, OPTIONS, OPTS)

%%MODEL: amin TV(S) + ||S||_1 + ||E||_1
%% Global constants and defaults

if isfield(OPTS,'RHO'),
    rho = OPTS.RHO;
else
    rho = 1.8;
end

if isfield(OPTS,'TOL'),
    TOL = OPTS.TOL;
else
    TOL = 1e-7;
end

if isfield(OPTS,'MAX_ITER'),
    MAX_ITER = OPTS.MAX_ITER;
else
    MAX_ITER = 50;
end

%% Set up ADMM

[m, n] = size(U);
[dim,dim2] = size(C);

w = zeros(n,1);
s = zeros(m,1);
xi = zeros(m,1);
e = zeros(m,1);
z = zeros(dim,1);

y1 = zeros(dim,1);
y2 = zeros(m,1);
y3 = zeros(m,1);

mu1 = 1.25/norm(v);
mu2 = mu1;
mu3 = mu1;

rho = 1.8;

% mu1 = 1e-4;
% mu2 = mu1;
% mu3 = mu2;

% precompute static variables for a-update (projection on to Ua=v-s)
P = (U'*U) \ (U');

fftUpdate = @(x) real(vec(ifft2(fft2(reshape(x,[M, N, channels])) .* Bf)));
C0 = C(:,idx);

%% ADMM solver

converge = false;
iter     = 0;

weights1 = ones(dim2,1);
weights1(1:M-1:end) = 0;
weights1(1) = 1;

weights2 = ones(dim2,1);
weights2(end-N:end) = 0;
weights = [weights1;weights2];
% weights = ones(dim,1);

B = zeros(dim2,1);
while ~converge && iter < MAX_ITER
    iter = iter + 1;
    
    % w update
    w = P * (v - xi - e - y2/mu2);
    
    Uw = U*w;
    r = v - Uw - y3/mu3;
    
    %xi update
    xi = shrinkage(0.5*(r - e + s - y2/mu2),lambdaS/(mu2+mu3));
    
    %e update
    e = shrinkage(r - xi, lambdaE/mu3);
    
    %s update
    B(idx) = (xi + y2/mu2);
    s0 = fftUpdate(C'*(z - y1/mu1) + B);
    s = s0(idx);
    
    %z update
    Cs = C0*s;
    z = shrinkage(Cs + y1/mu1, lambdaZ/mu1*weights);

    % dual1 update
    h1 = Cs - z;
    y1 = y1 + mu1 * h1;
    
    % dual2 update
    h2 = xi - s;
    y2 = y2 + mu2 * h2;
    
    % dual3 update
    h3 = Uw + e + xi - v;
    y3 = y3 + mu2 * h3;
    
    mu1 = min(rho * mu1, 1e10); 
    mu2 = min(rho * mu2, 1e10);
    mu3 = min(rho * mu3, 1e10);
        
    % diagnostics, reporting, termination checks
    
    if (norm(h1) < TOL && norm(h2) < TOL && norm(h3) < TOL)
        converge = true;    
    end
   
end

end
%%

function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
end



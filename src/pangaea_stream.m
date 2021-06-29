%%Kyle Gilman July 2019
%%Some code adopted from Jun He, Laura Balzano, and Arthur Szlam 2012

function [ Unew, STATUSnew, OPTSnew ] = pangaea_stream( y_Omega, idx, mask, U0, Bf, C, STATUS, OPTIONS, OPTS)


if isfield(OPTIONS,'DIM_M'),
    DIM_M = OPTIONS.DIM_M;
else
    error('Should specify OPTIONS.DIM_M data ambient dimension!!!\n');
end

if isfield(OPTIONS,'ITER_MIN'),
    MIN_ITER = OPTIONS.ITER_MIN;
else
    MIN_ITER = 5;
end

if isfield(OPTIONS,'ITER_MAX'),
    ITER_MAX = OPTIONS.ITER_MAX;
else
    ITER_MAX = 60;
end

if isfield(OPTIONS,'TOL'),
    TOL = OPTIONS.TOL;
else
    TOL = 1e-6;
end


if isfield(OPTIONS,'RANK'),
    RANK = OPTIONS.RANK;
else
    error('Should specify OPTIONS.RANK!!!\n');
end

% initilization 
if STATUS.init == 0,    
    STATUS.init         = 1; % Do not enter this initial part any longer   
    STATUS.curr_iter    = 0; % For debug

    OPTS.TOL        = TOL;
    OPTS.MAX_ITER   = MIN_ITER;  % the max iteration of ADMM at level=0   

    if isfield(OPTIONS,'rho'),
        OPTS.RHO = OPTIONS.rho;
    else
        OPTS.RHO = 1.8;
    end
    
    U0 = orth(randn(DIM_M,RANK));
end

%%%%%%%%%%%%%%%
% main framework of GRASTA
U_Omega = U0(idx,:);

% % ADMM OF PARAMETERS

[e_t, s_t, xi_t, w,ldual3] = admm_srp_pangaea(U_Omega, y_Omega, Bf, C, idx, OPTIONS, OPTS);  %MODEL: TV(S) + ||S||_1 + ||E||_1

gamma_1 = ldual3 + OPTS.RHO*(U_Omega*w + e_t + xi_t - y_Omega);

UtDual_omega = U_Omega' * gamma_1;
gamma_2 = U0 * UtDual_omega;
gamma = zeros(DIM_M,1);
gamma(idx) = gamma_1;
gamma = gamma - gamma_2;

gamma_norm = norm(gamma);
w_norm     = norm(w);

% 4. update the gradient for further step size update
STATUS.last_gamma  = gamma;
STATUS.last_w      = w;

% Take the gradient step along Grassmannian geodesic.
alpha = w/w_norm;
beta  = gamma/gamma_norm;

sG = norm(gamma)*norm(w);

t = OPTIONS.step / sG;

step = (cos(sG*t)-1)*U0*(alpha*alpha')  - sin(sG*t)*beta*alpha';

U0 = U0 + step;

%%

STATUS.e_t   = e_t;
STATUS.s_t   = s_t;
STATUS.xi_t  = xi_t;
STATUS.w     = w;
STATUS.SCALE = 1;
STATUS.curr_iter = STATUS.curr_iter + 1;

STATUS.step = t;
%%%%%%%%%%%%%%%%%%%%%%

Unew = U0;
STATUSnew = STATUS;
OPTSnew = OPTS;
end


function [Ltil, Stil] = run_grasta(Y, varargin)

[Mask, height, width, max_cycles, startStep, r, admm_iters] = parseInputs(Y,varargin{:});

[dim,frame_count] = size(Y);
U = orth(randn(dim,r));
M = height;
N = width;

OPTS            = struct(); % initial a empty struct for OPTS
status.init     = 0;   % status of grasta at each iteration

OPTIONS.DIM_M   = dim; % video ambient dimension

% GRASTA parameters
OPTIONS.RANK                = r;  % the estimated low-rank
OPTIONS.rho                 = 1.8;    
OPTIONS.ITER_MAX            = 20; 
OPTIONS.ITER_MIN            = admm_iters;    % the min iteration allowed for ADMM at the beginning

OPTIONS.USE_MEX             = 0;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.                                     

% Initialize the rough subspace
OPTIONS.CONSTANT_STEP       = 0.1;   % use adaptive step-size to initialize the subspace
OPTIONS.MAX_LEVEL           = 20;
OPTIONS.MAX_MU              = 10000; % set max_mu large enough for initial subspace training
OPTIONS.MIN_MU              = 1;


%Step and diagnostics initialization
t = startStep;
obj_log = zeros(frame_count*max_cycles,1);

Stil = zeros(dim,frame_count);
Ltil = zeros(dim,frame_count);

%Begin stochastic gradient descent
iter = 1;
for outiter = 1 : max_cycles
    
    t = t /outiter;
    OPTIONS.CONSTANT_STEP = t;
    tStart = tic;
    
%     p = randperm(frame_count);
    p = 1:frame_count;
    
    for i=1:frame_count
%        idx = find(Mask(:,p(i))==1);
        rho = 0.2;
        mask_idx = find(Mask(:,p(i)));
        idx = randsample(mask_idx,round(rho*length(mask_idx)));
       
       y = Y(:,p(i));
       y_Omega = y(idx);
%        y_Omega = y_Omega / norm(y_Omega);

       [U, status, OPTS] = grasta_stream(y_Omega, idx, U, status, OPTIONS, OPTS);
        
       obj_log(iter) = norm(y - U * status.w * status.SCALE, 1);
       
       Ltil(:,p(i)) = U*status.w;
       
       Stil(:,p(i)) = y - U * status.w * status.SCALE;
     
       iter = iter + 1;
    end
    
    
    t_end = toc(tStart);
    
    fprintf('Training %d/%d: %.2f seconds, %.2f fps, grasta_t %.2e \n',...
        outiter, max_cycles,t_end, frame_count/t_end,status.grasta_t);

end

obj_log = obj_log(1:iter-1);
% figure,semilogy(obj_log); title('TV-GRASTA Objective Fxn'); xlabel('Iteration')

end



function [M, height, width, max_cycles, startStep, rank, admm_iters] = parseInputs(Y,opts)
if ~exist('opts','var')
    opts = struct();
end
max_cycles          = parseField(opts,'max_cycles',5);
M                   = parseField(opts,'M',1);
startStep           = parseField(opts,'startStep',0.5);
rank                = parseField(opts,'rank',1);
admm_iters          = parseField(opts,'admm_iters',50);
height              = parseField(opts,'height');
width               = parseField(opts,'width');

end

function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end
end
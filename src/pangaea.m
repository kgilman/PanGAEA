function [Ltil, Etil, Stil] = pangaea(Y,lambdaZ,lambdaE, lambdaS, varargin)

[Mask, height, width, channels, max_cycles, startStep, r, admm_iters] = parseInputs(Y,varargin{:});

[dim,frame_count] = size(Y);
U = orth(randn(dim,r));
M = height;
N = width;

OPTS            = struct(); % initial a empty struct for OPTS
status.init     = 0;   % status at each iteration

OPTIONS.DIM_M   = dim; % video ambient dimension

% PanGAEA parameters
OPTIONS.RANK                = r;  % the estimated low-rank
OPTIONS.rho                 = 1.8;    
OPTIONS.ITER_MAX            = 100; 
OPTIONS.ITER_MIN            = admm_iters;    % the min iteration allowed for ADMM at the beginning   
OPTIONS.lambdaZ             = lambdaZ;
OPTIONS.lambdaE             = lambdaE;
OPTIONS.lambdaS             = lambdaS;
OPTIONS.M                   = M;
OPTIONS.N                   = N;
OPTIONS.channels            = channels;

Etil = zeros(dim,frame_count);
Stil = zeros(dim,frame_count);
Ltil = zeros(dim,frame_count);

%TV Norm setup
[C2, ~] = Cc(M,N,1,2);
Hf = fft2(reshape(full(C2' * C2(:,1)),[M, N]));
Bf = repmat(1 ./ (1 + Hf),1,1,channels);

[C,~] = Cc(M,N,channels,2);

t = startStep;  %starting step size

%Begin stochastic gradient descent
iter = 1;
for outiter = 1 : max_cycles
    
    times = zeros(frame_count,1);
    
    %%Step-size attenuation policy
    t = t /outiter;
%     t = t /10;
    OPTIONS.step = t;
   
    %%Frame-order policy
    p = randperm(frame_count);
%     p = 1:frame_count;
    
    for i=1:frame_count
       
       if(outiter < max_cycles)
           rho = 0.2;
           mask_idx = find(Mask(:,p(i)));
           idx = randsample(mask_idx,round(rho*length(mask_idx)));
       else
           idx = find(Mask(:,p(i)));
       end
%        idx = find(Mask(:,p(i)));
       y = Y(:,p(i));
       y_Omega = y(idx);

       tStart = tic;
       [U, status, OPTS] = pangaea_stream(y_Omega, idx, Mask(:,p(i)), U, Bf, C, status, OPTIONS, OPTS);
       times(iter) = toc(tStart);
       
       % Store estimated variables
       Ltil(:,p(i)) = U*status.w;   %low-rank background
       
       e0 = zeros(dim,1);
       e0(idx) = status.e_t;
       Etil(:,p(i)) = e0;           % sparse corruptions
       
       s0 = zeros(dim,1);
       s0(idx) = status.xi_t;
       Stil(:,p(i)) = s0;           % sparse fg
     
       iter = iter + 1;
    end
    
    tTotal = sum(times);
    fprintf('Training %d/%d: %.2f seconds, %.2f fps, step %.2e \n',...
        outiter, max_cycles,tTotal, frame_count/tTotal,status.step);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Internal Functions%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circulant 1D first differences matrix
function D = Dc(n)
D = spdiags([-ones(n,1), ones(n,1)],[0, 1],n,n);
D(n,1) = 1;
end

function [C, osiz] = Cc(m,n,p,dim)
m = double(m);
n = double(n);
C = [kron(speye(p),kron(speye(n),Dc(m)));         % columns
     kron(speye(p),kron(Dc(n),speye(m)))];        % rows
if dim == 3
    C = [C; kron(Dc(p),kron(speye(n),speye(m)))]; % frames
end
osiz = [size(C,1), 1];

end


function [M, height, width, channels, max_cycles, startStep, rank, admm_iters] = parseInputs(Y,opts)
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
channels               = parseField(opts,'channels');

end

function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

end
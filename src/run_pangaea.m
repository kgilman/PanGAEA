function [pano, L, E, S, S_disp, Lreg, Ereg, Sreg_disp] = run_pangaea(Yreg,Ynoisy,varargin)
%              
% Description:  Performs foreground-background separation via the
%               PanGAEA method
%
% Inputs:       Yreg: registered, vectorized panoramic frames M*N x
%               num_frames
%               Ynoisy: original frames array of size m x n x num_frames
%
% Kyle Gilman, University of Michigan
% July 2019
%
% Some code adopted from B. E. Moore, C. Gao, and R. R. Nadakuditi


% Parse inputs
[mask, height, width, lambdaZ, lambdaE, lambdaS, max_cycles, startStep, rank, admm_iters, cleanLS] = parseRPCAInputs(Yreg,varargin{:});
isRGB = (size(Ynoisy,3) > 1);

% Perform TV-GRASTA
opts.M = mask;
opts.max_cycles = max_cycles;
opts.startStep = startStep;
opts.rank = rank;
opts.admm_iters = admm_iters;
opts.height = height;
opts.width = width;
if(isRGB)
    opts.channels = 3;
else
    opts.channels = 1;
end
[Lhat, Ehat, Shat] = pangaea(Yreg,lambdaZ,lambdaE,lambdaS, opts);

if isRGB
    %
    % Color images
    %
    
    % Generate forground/background video
    Lreg = reshape(Lhat,[height, width, 3, size(Ynoisy,4)]);
    Ereg = reshape(Ehat,[height, width, 3, size(Ynoisy,4)]);
    Sreg = reshape(Shat,[height, width, 3, size(Ynoisy,4)]);
    M = logical(reshape(mask,[height, width, 3, size(Ynoisy,4)]));

    if cleanLS
        [Lreg, Sreg] = adjustLS(Lreg,Sreg,M);
    end
    
    [Lreg, Ereg, Sreg_disp] = formatForDisplay(Lreg,Ereg,Sreg,M);
    
    L = pano2video_RGB(Lreg,mask,height,width,size(Ynoisy));
    E = pano2video_RGB(Ereg,mask,height,width,size(Ynoisy));
    S = pano2video_RGB(Sreg,mask,height,width,size(Ynoisy));
    S_disp = pano2video_RGB(Sreg_disp,mask,height,width,size(Ynoisy));
    
    % Generate panorama
    [uhat, ~, ~] = svds(reshape(Lreg,[],size(Ynoisy,4)),1);
    pano = reshape(uhat,height,width,3);
    pano = uint8(round(pano * 255 / max(pano(:))));
else
    %
    % Grayscale images
    %
    
    % Generate forground/background video
    Lreg = reshape(Lhat,[height, width, size(Yreg,2)]);
    Ereg = reshape(Ehat,[height, width, size(Yreg,2)]);
    Sreg = reshape(Shat,[height, width, size(Yreg,2)]);
    M = logical(reshape(mask,[height, width, size(Yreg,2)]));

    if cleanLS
        [Lreg, Sreg] = adjustLS(Lreg,Sreg,M);
        Lreg = cleanBackground(Lreg,M);
    end

    [Lreg, Ereg, Sreg_disp] = formatForDisplay(Lreg,Ereg,Sreg,M);
    L = pano2video(Lreg,mask,height,width,size(Ynoisy));
    E = pano2video(Ereg,mask,height,width,size(Ynoisy));
    S = pano2video(Sreg,mask,height,width,size(Ynoisy));
    S_disp = pano2video(Sreg_disp,mask,height,width,size(Ynoisy));
    % Generate panorama
    [uhat, ~, ~] = svds(reshape(Lreg,[],size(Ynoisy,4)),1);
    pano = reshape(uhat,[height, width]);
end

% Parse registration inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [method, T] = parseRegistrationInputs(opts)
if ~exist('opts','var')
    opts = struct();
end
method   = parseField(opts,'method','temporal');
T        = parseField(opts,'T',nan);


% Parse RPCA inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask, height, width, lambdaZ, lambdaE, lambdaS, max_cycles, startStep, rank, admm_iters, cleanLS] = parseRPCAInputs(Yreg,opts)
if ~exist('opts','var')
    opts = struct();
end
mask            = parseField(opts,'mask');
height          = parseField(opts,'height');
width           = parseField(opts,'width');
lambdaZ         = parseField(opts,'lambdaZ',1);
lambdaE         = parseField(opts,'lambdaE',1.5);
lambdaS         = parseField(opts,'lambdaS',1.5);
max_cycles      = parseField(opts,'max_cycles',5);
startStep       = parseField(opts,'startStep',0.5);
rank            = parseField(opts,'rank',1);
admm_iters      = parseField(opts,'admm_iters',100);
cleanLS         = parseField(opts,'cleanLS',true);


% Parse struct field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

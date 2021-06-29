function [pano, L, S, S_disp, Lreg, Sreg, Sreg_disp] = GRASTA(Yreg,Ynoisy,varargin)


% Parse RPCA inputs
[mask, height, width, max_cycles, startStep, rank, admm_iters, cleanLS] = parseRPCAInputs(Yreg,varargin{:});
% isRGB = ndims(Yreg) == 4;
isRGB = (size(Ynoisy,3) > 1);

% Perform GRASTA
opts.M = mask;
opts.max_cycles = max_cycles;
opts.startStep = startStep;
opts.rank = rank;
opts.admm_iters = admm_iters;
opts.height = height;
opts.width = width;
[Lhat, Shat] = run_grasta(Yreg, opts);

if isRGB
    %
    % Color images
    %
    
    % Generate forground/background video
    Lreg = reshape(Lhat,[height, width, 3, size(Yreg,2)]);
    Sreg = reshape(Shat,[height, width, 3, size(Yreg,2)]);
    M = logical(reshape(mask,[height, width, 3, size(Yreg,2)]));

    if cleanLS
        [Lreg, Sreg] = adjustLS(Lreg,Sreg,M);
    end
    
    [Lreg, ~, Sreg_disp] = formatForDisplay(Lreg,[],Sreg,M);
    
    L = pano2video_RGB(Lreg,mask,height,width,size(Ynoisy));
    S = pano2video_RGB(Sreg,mask,height,width,size(Ynoisy));
    S_disp = pano2video_RGB(Sreg_disp,mask,height,width,size(Ynoisy));

    % Generate panorama
    [uhat, ~, ~] = svds(reshape(Lreg,[],size(Yreg,2)),1);
    pano = reshape(uhat,height,width,3);
    pano = uint8(round(pano * 255 / max(pano(:))));
else
    %
    % Grayscale images
    %
    
    % Generate forground/background video
    Lreg = reshape(Lhat,[height, width, size(Yreg,2)]);
    Sreg = reshape(Shat,[height, width, size(Yreg,2)]);
    M = logical(reshape(mask,[height, width, size(Yreg,2)]));

    if cleanLS
        [Lreg, Sreg] = adjustLS(Lreg,Sreg,M);
        Lreg = cleanBackground(Lreg,M);
    end

    [Lreg, ~, Sreg_disp] = formatForDisplay(Lreg,[],Sreg,M);
    L = pano2video(Lreg,mask,height,width,size(Ynoisy));
    S = pano2video(Sreg,mask,height,width,size(Ynoisy));
    S_disp = pano2video(Sreg_disp,mask,height,width,size(Ynoisy));

    % Generate panorama
    [uhat, ~, ~] = svds(reshape(Lreg,[],size(Ynoisy,4)),1);
    pano = reshape(uhat,[height, width]);
end


end

% Parse RPCA inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask, height, width, max_cycles, startStep, rank, admm_iters, cleanLS] = parseRPCAInputs(Yreg,opts)
if ~exist('opts','var')
    opts = struct();
end
mask            = parseField(opts,'mask');
height          = parseField(opts,'height');
width           = parseField(opts,'width');
max_cycles      = parseField(opts,'max_cycles',5);
startStep       = parseField(opts,'startStep',0.5);
rank            = parseField(opts,'rank',1);
admm_iters      = parseField(opts,'admm_iters',50);
cleanLS         = parseField(opts,'cleanLS',true);
end

% Parse struct field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

end
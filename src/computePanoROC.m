function roc = computePanoROC(E,VIDEOPATH,isRGB,scale)

%E is the mxnxp array of foreground images from PRPCA algorithms
%VIDEOPATH is the path to the ground truth binary images in DAVIS
%scale is the scaling of the videos with which to downsample the frames

filename = VIDEOPATH;
imagefiles = dir(filename);
imagefiles(1:2) = [];
nfiles = length(imagefiles);    % Number of files found
    
for ii=1:nfiles
    currentfilename = imagefiles(ii).name;
    currentimage = imread(strcat(filename,currentfilename));
    currentimage = imresize(currentimage,scale);
    binary = zeros(size(currentimage));
    binary(currentimage>0) = 1;
    if(ii==1)
        [m,n] = size(binary);
        gndTruth = zeros(m,n,nfiles);
    end
    gndTruth(:,:,ii) = binary;
end


gammas = linspace(0,1,100);

fp_tot = zeros(1,length(gammas));
tp_tot = zeros(1,length(gammas));
roc = zeros(length(gammas),2);
Z = 0;
tidx_total = 0;

for frame=1:size(E,4)
    bin_im = gndTruth(:,:,frame);
    truth_idx = find(bin_im);
    Z = Z + m*n;
    tidx_total = tidx_total + length(truth_idx);
    
    
    if(isRGB)
        Ef = sqrt(E(:,:,1,frame).^2 + E(:,:,2,frame).^2 + E(:,:,3,frame).^2);
    else
        Ef = E(:,:,1,frame);
    end
    
    for t=1:length(gammas)
        
        gam = gammas(t);
        Ethresh = fg_thresholding(Ef,gam);
        detected_idx = find(Ethresh);
        
        false_pos = setdiff(detected_idx,truth_idx);
        true_pos = intersect(truth_idx,detected_idx);
        
        fp_tot(t) = fp_tot(t) + length(false_pos);
        tp_tot(t) = tp_tot(t) + length(true_pos);
    end

    
end

for t=1:length(gammas)
    roc(t,1) = fp_tot(t)/(Z-tidx_total);
    roc(t,2) = tp_tot(t)/tidx_total;
end

end
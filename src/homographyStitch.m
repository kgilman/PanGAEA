function [Y,Ytrue,Ynoisy,Mask,M,N,STITCHING] = homographyStitch(VIDEOPATH,isRGB,scale,noise,noise_level)

filename = VIDEOPATH;
imagefiles = dir(filename);
nfiles = length(imagefiles);    % Number of files found
%nfiles = round(nfiles/2);

%GET THE ANCHOR FRAME
if(mod(nfiles,2)==0); nfiles = nfiles - 1; end
anchor_idx = (nfiles - 1)/2;

SURFthresh = 10;
%PERFORM HOMOGRAPHY REGISTRATION
frameidx = 1;
for ii=3:nfiles-1
    
    currentfilename = imagefiles(ii).name;
    currentimage = imread(strcat(filename,currentfilename));
    frameGray = single(rgb2gray(currentimage));
    frameGray = imresize(frameGray,scale);
    
    nextfilename = imagefiles(ii+1).name;
    nextimage = imread(strcat(filename,nextfilename));
    frameGray_next = single(rgb2gray(nextimage));
    frameGray_next = imresize(frameGray_next,scale);
    
    points1 = detectSURFFeatures(uint8(frameGray),'MetricThreshold',SURFthresh);
    points2 = detectSURFFeatures(uint8(frameGray_next),'MetricThreshold',SURFthresh);
    
    [features1,vpts1] = extractFeatures(uint8(frameGray),points1);
    [features2,vpts2] = extractFeatures(uint8(frameGray_next),points2);
    
    indexPairs = matchFeatures(features1,features2,'Unique',true);

    matchedPoints1 = vpts1(indexPairs(:,1));
    matchedPoints2 = vpts2(indexPairs(:,2));
    
    
    tform = estimateGeometricTransform(matchedPoints1,matchedPoints2,'projective','Confidence',99.9,'MaxNumTrials',2000);
    
    H_MAT{frameidx} = tform.T;
    frameidx = frameidx + 1;
    
    if(isRGB)
        store_frame = single(imresize(currentimage,scale));
    else
        store_frame = frameGray;
    end
    
    Ytrue(:,:,:,ii-2) = store_frame;
    
    if(noise)
        store_frame = single(imnoise(uint8(store_frame),'salt & pepper',noise_level));
    end
    store_frame = double(store_frame) / 255;

    images{ii-2} = store_frame;
    Ynoisy(:,:,:,ii-2) = store_frame;
    
end

if(isRGB)
    store_frame = single(imresize(nextimage,scale));
else
    store_frame = frameGray_next;
end

if(noise)
    store_frame = single(imnoise(uint8(store_frame),'salt & pepper',noise_level));
end

store_frame = double(store_frame) / 255;
    
images{ii-1} = store_frame;
Ynoisy(:,:,:,ii-1) = store_frame;


% %PLAY THE NORMAL VIDEO
% for i=1:frameidx-1
%     figure(1),imagesc(images{i}); axis image
% end


%COMPUTE GLOBAL HOMOGRAPHIES FOR K<K~
for frame=1:anchor_idx-1
    
    H = eye(3,3);
    for hidx = frame:anchor_idx-1
        H = H_MAT{hidx}*H;
    end
    
    H_global{frame} = H;
 
end

%COMPUTE GLOBAL HOMOGRAPHIES FOR K>=K~
for frame=anchor_idx:frameidx
    
    H = eye(3,3);
    for hidx = anchor_idx:frame-1
        H = H*H_MAT{hidx}';
    end
    
    H_global{frame} = inv(H');
    
end

%CONSTRUCT THE STITCHING STRUCTURE
for i=1:length(H_global)
    It(i).image = images{i};
    It(i).tform = H_global{i};
end

%COMPUTE STITCHED VIDEO
%[Is,frames,alpha] = stitchImages(It,'mode','absolute');
[Is, frames, Ias, coord, indices, alpha] = stitchImages2(It,'mode','absolute');
    
%SHOW STITCHED PANORAMA
% figure,imagesc(Is); axis image

% %SHOW BORDERS
% figure(5),imagesc(alpha); axis image

[M,N,c] = size(Is);

%PLAY SITCHED VIDEO, PREPARE THE FRAMES AND CONSTRUCT OBSERVATION MATRIX
t = length(H_global);

Y = zeros(M*N*c,t);
for i=1:t
    im = frames{i};
    figure(6),imagesc((im));  axis image
    Y(:,i) = im(:);

end


[m_o,n_o,c] = size(store_frame);

Mask = ~isnan(Y);
Y(isnan(Y)) = 0;

STITCHING.Ias = Ias;
STITCHING.coord = coord;
STITCHING.indices = indices;
STITCHING.It = It;
STITCHING.anchor = anchor_idx;
STITCHING.m_o = m_o;
STITCHING.n_o = n_o;

end

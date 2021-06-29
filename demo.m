%% Demo script of PanGAEA, GRASTA, RPCA, and PRPCA on video foreground-background separation of DAVIS Challenge videos
%
% Kyle Gilman, July 2019
%
% RPCA and PRPCA code and some of top-level scripts from B. E. Moore, C.
% Gao, and R. R. Nadakuditi 2017
%
% GRASTA code from Jun He, Laura Balzano, and Arthur Szlam 2012
%

close all
clear
clc
addpath('src/');

% Load Panoramic Video Data

rng(1);


%$SET THE VIDEO PATH - change this to point to the dataset!

video_path = '/Users/kgilman/Desktop/Desktop - kgilman-mbp/Datasets/DAVIS-3/JPEGImages/480p/';


%%Choose a video

% dataset = 'dog-gooses/';
% dataset = 'horsejump-high';
% dataset = 'horsejump-low';
% dataset = 'lucia';
dataset = 'paragliding';
% dataset = 'swing';
% dataset = 'tennis';
% dataset = 'flamingo';
% dataset = 'paragliding-launch';
% dataset = 'stroller';
% dataset = 'car-roundabout';
% dataset = 'car-shadow';
% dataset = 'hockey';
% dataset = 'blackswan';
% dataset = 'dance-jump';
% dataset = 'hike';
% dataset = 'bmx-trees';
% dataset = 'dance-twirl';

VIDEOPATH = strcat(video_path,dataset,'/');

%Also change below path!
video_path_truth = '/Users/kgilman/Desktop/Desktop - kgilman-mbp/Datasets/DAVIS-3/Annotations/480p/';
VIDEOPATH_TRUTH = strcat(video_path_truth,dataset,'/');

scale = 0.25;   %resolution
isRGB = 1;          %1 for RGB; 0 for grayscale
noise = 0;      %1 for impulse noise; 0 for clean video
noise_level = 0.2;  %Bernoulli-p of impulse noise

[Yreg,Ytrue,Ynoisy,mask,height,width,~] = homographyStitch(VIDEOPATH,isRGB,scale,noise,noise_level);

% Normalize the data to be in [0 1] range
Ytrue = Ytrue - min(min(min(min(Ytrue))));
Ytrue = Ytrue / max(max(max(max(Ytrue))));

Ynoisy = Ynoisy - min(min(min(min(Ynoisy))));
Ynoisy = Ynoisy / max(max(max(max(Ynoisy))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Noisy (PanGAEA 1)

opts = struct();

opts.lambdaZ = 1;
opts.lambdaS = 0.5;
opts.lambdaE = 1;

opts.cleanLS = false;
opts.max_cycles = 7;

opts.mask = mask;
opts.height = height;
opts.width = width;
opts.startStep = 1;
tStart = tic;
[pano, L_tvgrasta, E_tvgrasta, S_tvgrasta, S_tvgrasta_disp, Lreg_tvgrasta, Ereg_tvgrasta, Sreg_tvgrasta] = run_pangaea(Yreg,Ynoisy,opts);
tElasped_tv_grasta = toc(tStart);

roc_tvgrasta = computePanoROC(S_tvgrasta,VIDEOPATH_TRUTH,isRGB,scale);
area_tvgrasta = abs(trapz(roc_tvgrasta(:,1), roc_tvgrasta(:,2)));
figure, plot(roc_tvgrasta(:,1),roc_tvgrasta(:,2));

%%
%%Visualize registered decomposition
Yreg_disp = reshape(Yreg,size(Lreg_tvgrasta));
% movie = [Yreg_disp; Lreg_tvgrasta; Sreg_tvgrasta; Ereg_tvgrasta];
movie = [Yreg_disp; Lreg_tvgrasta; Sreg_tvgrasta];
opts2 = struct();
% opts2.ylabels = {'Sparse', 'Foreground', 'Background','Observed'};
PlayMovie(movie,opts2);

% movie = [squeeze(Ynoisy); squeeze(L_tvgrasta); squeeze(S_tvgrasta_disp); squeeze(E_tvgrasta)];
movie = [squeeze(Ynoisy); squeeze(L_tvgrasta); squeeze(S_tvgrasta_disp)];
opts3 = struct();
% opts3.ylabels = {'Sparse', 'Foreground', 'Background', 'Observed'};
opts3.ylabels = {'Foreground', 'Background', 'Observed'};
PlayMovie(movie,opts3);

% Visualize panorama
figure();
imshow(pano,[]);
title('Panoramic background: TV-GRASTA');

% Visualize key frames of decomposition
LS = L_tvgrasta + S_tvgrasta_disp;
% idx1 = 1;
% idx2 = 17;
% idx1 = 10;
% idx2 = 20;
% idx3 = 30;
% idx4 = 40;

idx1 = 15;
idx2 = 30;
idx3 = 45;
idx4 = 60;

% keyframes = [
%     [Ytrue(:,:,:,idx1);  Ynoisy(:,:,:,idx1);  L_tvgrasta(:,:,:,idx1) ; E_tvgrasta(:,:,:,idx1); S_tvgrasta_disp(:,:,:,idx1)] ...
%     [Ytrue(:,:,:,idx2); Ynoisy(:,:,:,idx2); L_tvgrasta(:,:,:,idx2); E_tvgrasta(:,:,:,idx2); S_tvgrasta_disp(:,:,:,idx2)]
% ];

keyframes = [
    [Ytrue(:,:,:,idx1);   L_tvgrasta(:,:,:,idx1) ; S_tvgrasta_disp(:,:,:,idx1)] ...
    [Ytrue(:,:,:,idx2);  L_tvgrasta(:,:,:,idx2); S_tvgrasta_disp(:,:,:,idx2)] ...
    [Ytrue(:,:,:,idx3);  L_tvgrasta(:,:,:,idx3); S_tvgrasta_disp(:,:,:,idx3)]...
    [Ytrue(:,:,:,idx4);  L_tvgrasta(:,:,:,idx4); S_tvgrasta_disp(:,:,:,idx4)]
];

figure();
imagesc(keyframes); axis image;  colormap gray; caxis([-1 1])
title('Key frames of decomposition');
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GRASTA

opts = struct();
opts.cleanLS = false;
opts.max_cycles = 10;
           
opts.mask = mask;
opts.height = height;
opts.width = width;
opts.startStep = 1;
tStart = tic;
% [pano, L_grasta, S_grasta, Lreg_grasta, Sreg_grasta] = GRASTA(Yreg,Ynoisy,opts);
[pano, L_grasta, S_grasta, S_grasta_disp, Lreg_grasta, Sreg_grasta, Sreg_grasta_disp] = GRASTA(Yreg,Ynoisy,opts);
tElapsed_grasta = toc(tStart);
roc_grasta = computePanoROC(S_grasta,VIDEOPATH_TRUTH,isRGB,scale);
area_grasta = abs(trapz(roc_grasta(:,1), roc_grasta(:,2)));
figure,plot(roc_grasta(:,1),roc_grasta(:,2));
title('GRASTA')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%Visualize registered decomposition 
Yreg_disp = reshape(Yreg,size(Lreg_grasta));
movie = [Yreg_disp; Lreg_grasta; Sreg_grasta];
opts2 = struct();
% opts2.ylabels = {'Sparse', 'Foreground', 'Background','Observed'};
PlayMovie(movie,opts2);

movie = [squeeze(Ynoisy); squeeze(L_grasta); squeeze(S_grasta_disp)];
opts3 = struct();
opts3.ylabels = {'Sparse', 'Foreground', 'Background', 'Observed'};
PlayMovie(movie,opts3);

% Visualize panorama
figure();
imshow(pano,[]);
title('Panoramic background: TV-GRASTA');

% Visualize key frames of decomposition
LS = L_grasta + S_grasta_disp;
% idx1 = 1;
% idx2 = 17;
idx1 = 10;
idx2 = 30;

keyframes = [
    [Ytrue(:,:,:,idx1);  Ynoisy(:,:,:,idx1);  L_grasta(:,:,:,idx1); S_grasta_disp(:,:,:,idx1)] ...
    [Ytrue(:,:,:,idx2); Ynoisy(:,:,:,idx2); L_grasta(:,:,:,idx2); S_grasta_disp(:,:,:,idx2)]
];

figure();
imagesc(keyframes); axis image;  colormap gray; caxis([-1 1])
title('Key frames of decomposition');
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% Noiseless RPCA

% Perform Noiseless PRPCA
opts = struct();
opts.nIters = 25;
opts.mask = mask;
opts.height = height;
opts.width = width;
opts.cleanLS = false;
tStart = tic;
[pano_rpca, L_rpca, S_rpca, Lreg_rpca, Sreg_rpca] = PRPCA_noiseless(Yreg,Ynoisy,opts);
tElapsed_rpca = toc(tStart);
%Compute ROC Curve
roc_prpca_noiseless = computePanoROC(S_rpca,VIDEOPATH_TRUTH,isRGB,scale);
area = abs(trapz(roc_prpca_noiseless(:,1), roc_prpca_noiseless(:,2)));
figure,plot(roc_prpca_noiseless(:,1),roc_prpca_noiseless(:,2));
title('RPCA')
% 
%
%% PRPCA

rng(1);

% Perform PRPCA
opts = struct();
opts.mask = mask;
opts.height = height;
opts.width = width;
opts.nIters = 150;
opts.nItersS = 10;
opts.cleanLS = false;
opts.lambdaS = 0.02;
opts.lambdaE = 0.02;
tStart = tic;
[pano_prpca, L_prpca, E_prpca, S_prpca, S_prpca_disp, Lreg_prpca, Ereg_prpca, Sreg_prpca, Sreg_prpca_disp] = PRPCA(Yreg,Ynoisy,opts);
tElapsed_prpca = toc(tStart);

% load(strcat('/Users/kgilman/Desktop/Desktop - kgilman-mbp/Quals/panoramicRPCA-master/prpca_results_',dataset,'_rgb.mat'));

%% 
% %Compute ROC Curve
roc_prpca= computePanoROC(S_prpca,VIDEOPATH_TRUTH,isRGB,scale);
area = abs(trapz(roc_prpca(:,1), roc_prpca(:,2)));
figure,plot(roc_prpca(:,1),roc_prpca(:,2));
title('PRPCA')

% Visualize registered decomposition
movie = [Yreg_disp;Lreg_prpca; Sreg_prpca_disp; Ereg_prpca];
opts2 = struct();
opts2.ylabels = {'sparse', 'foreground', 'background','observed'};
PlayMovie(movie,opts2);

% Visualize decomposition from original persepctive
movie = [Ynoisy; L_prpca; S_prpca_disp; E_prpca];
opts3 = struct();
opts3.ylabels = {'sparse', 'foreground', 'background', 'observations'};
PlayMovie(movie,opts3);

% Visualize panorama
figure();
imshow(pano_prpca,[]);
title('Panoramic background: PRPCA');

% Visualize key frames of decomposition
LS_prpca = L_prpca + S_prpca;

keyframes = [
    [Ytrue(:,:,:,idx1);  Ynoisy(:,:,:,idx1); L_prpca(:,:,:,idx1) ; E_prpca(:,:,:,idx1);  S_prpca_disp(:,:,:,idx1)] ...
    [Ytrue(:,:,:,idx2); Ynoisy(:,:,:,idx2); L_prpca(:,:,:,idx2); E_prpca(:,:,:,idx2); S_prpca_disp(:,:,:,idx2)]
];

figure();
imagesc(keyframes); axis image;  colormap gray; caxis([-1 1])
title('Key frames of decomposition');
% 
%% ROC CURVES

figure,
plot(roc_tvgrasta(:,1),roc_tvgrasta(:,2),':s','LineWidth',2);
hold on 
plot(roc_grasta(:,1),roc_grasta(:,2),':s','LineWidth',2);
hold on
plot(roc_prpca_noiseless(:,1),roc_prpca_noiseless(:,2),':s','LineWidth',2)
hold on;
plot(roc_prpca(:,1),roc_prpca(:,2),':s','LineWidth',2);

legend({'PanGAEA', 'GRASTA','RPCA','PRPCA'},'Location','SouthEast','FontSize',25);
xlabel('FP Rate')
ylabel('TP Rate')
title(strcat('DAVIS: ',dataset))
xticks(0:0.1:1);
yticks(0:0.1:1);
set(gca,'FontSize',15)

%% PSNR Measurements
filename_img = VIDEOPATH;
imagefiles_img = dir(filename_img);
imagefiles_img(1:2) = [];

filename_bin = VIDEOPATH_TRUTH;
imagefiles_bin = dir(filename_bin);
imagefiles_bin(1:2) = [];

nfiles = length(imagefiles_img);

psnr_prpca_log = zeros(nfiles,1);
psnr_tv_grasta_log = zeros(nfiles,1);
psnr_grasta_log = zeros(nfiles,1);
psnr_rpca_log = zeros(nfiles,1);

for i=1:size(Yreg,2)
%     img_rpca = adjustIm(L_rpca(:,:,:,i) + S_rpca(:,:,:,i));
%     img_prpca = adjustIm(L_prpca(:,:,:,i) + S_prpca(:,:,:,i));
    img_tv_grasta = adjustIm(L_tvgrasta(:,:,:,i) + S_tvgrasta(:,:,:,i));
%     img_grasta = adjustIm(L_grasta(:,:,:,i) + S_grasta(:,:,:,i));

    %Get the truth foreground image
    currentfilename = imagefiles_img(i).name;
    currentimage = imread(strcat(filename_img,currentfilename));
    currentimage = double(imresize(currentimage,scale));

    truth_img = currentimage;

%     psnr_rpca_log(i) = psnr(uint8(img_rpca),uint8(truth_img));
%     psnr_prpca_log(i) = psnr(uint8(img_prpca),uint8(truth_img));
    psnr_tv_grasta_log(i) = psnr(uint8(img_tv_grasta),uint8(truth_img));
%     psnr_grasta_log(i) = psnr(uint8(img_grasta),uint8(truth_img));

end

% mean_psnr_rpca = mean(psnr_rpca_log);
% mean_psnr_prpca = mean(psnr_prpca_log);
mean_psnr_tv_grasta = mean(psnr_tv_grasta_log);
% mean_psnr_grasta = mean(psnr_grasta_log);


idx = 40;

%Get the truth foreground image
currentfilename = imagefiles_img(idx).name;
currentimage = (imread(strcat(filename_img,currentfilename)));
currentimage = double(imresize(currentimage,scale));

truth_img = currentimage;

% img_rpca = adjustIm(L_rpca(:,:,:,idx) + S_rpca(:,:,:,idx));
% img_prpca = adjustIm(L_prpca(:,:,:,idx) + S_prpca(:,:,:,idx));
img_tv_grasta = adjustIm(L_tvgrasta(:,:,:,idx) + S_tvgrasta(:,:,:,idx));
% img_grasta = adjustIm(L_grasta(:,:,:,idx) + S_grasta(:,:,:,idx));

noisy_img = Ynoisy(:,:,:,idx);


figure,
subplot(4,2,1:2),imagesc(uint8(255*noisy_img)); axis image; title('Observed Image');
subplot(4,2,3),imagesc(uint8(truth_img)); axis image; title('True Foreground');
% subplot(4,2,4),imagesc(uint8(img_rpca)); axis image; title(strcat('RPCA PSNR: ',num2str(mean_psnr_rpca))); 
% subplot(4,2,5),imagesc(uint8(img_grasta));  axis image; title(strcat('GRASTA PSNR: ',num2str(mean_psnr_grasta))); 
% subplot(4,2,7),imagesc(uint8(img_prpca)); axis image; title(strcat('PRPCA PSNR: ',num2str(mean_psnr_prpca))); 
subplot(4,2,8),imagesc(uint8(img_tv_grasta)); axis image; title(strcat('TV-GRASTA PSNR: ',num2str(mean_psnr_tv_grasta))); 



function out = adjustIm(im)
    im = im - min(im(:));
    out = im / max(im(:)) * 255;
end


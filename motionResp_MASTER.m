%% Script to investigate the influence of respiration on motion correction
%
%   Written by Andrew S Bock Oct 2016

%% set defaults
figDir                          = '/Users/abock/MOTION_CORRECTION_TECHDEV/figures';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/heteroPhantom/042016/';
params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/heteroPhantom/051016';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/heteroPhantom/052016';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/homoPhantom/042016';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/HCLV1001/8792BT';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/FBIRN/052516';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/FBIRN/100416';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3001/081916a';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3002/082616a';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3003/090216';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3004/091916';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3005/092316';
params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/simulations';
params.despike                  = 0; % params.despike data
params.slicetiming              = 0; % do slice timing correction
params.refvol                   = 1; % reference volume = 1st TR
params.regFirst                 = 0; % register to the first run
runNum                          = 1;
b                               = find_bold(params.sessionDir);
Fs                              = 1/0.8; % Sampling frequency
%% Correct motion using mri_robust_register
motion_slice_correction(params,runNum);
%% Create reference volume (if referencing to first run)
% thisDir = fullfile(params.sessionDir,b{1});
% cd(thisDir);
% system('fslroi ./raw_f.nii.gz ./refVol.nii.gz 0 1');
% refFile = fullfile(thisDir,'refVol.nii.gz');
%% Correct motion using mcflirt
thisDir = fullfile(params.sessionDir,b{runNum});
cd(thisDir);
system(['mcflirt -in ./raw_f.nii.gz -out ./mcflirt.nii.gz -refvol 0 ' ...
    '-stats -mats -plots -report']);
%% Get motion params - mcflirt
clear mcflirtOut
thisDir = fullfile(params.sessionDir,b{runNum});
mats = listdir(fullfile(thisDir,'mcflirt.nii.gz.mat/MAT*'),'files');
cd(fullfile(thisDir,'mcflirt.nii.gz.mat'));
mcflirtOut = nan(length(mats),6);
for i = 1:length(mats)
    [x,y,z,pitch,yaw,roll] = convertMAT2tranrot(mats{i});
    mcflirtOut(i,1) = pitch*50;
    mcflirtOut(i,2) = yaw*50;
    mcflirtOut(i,3) = roll*50;
    mcflirtOut(i,4) = x;
    mcflirtOut(i,5) = y;
    mcflirtOut(i,6) = z;
end
%% Get motion params - robust register
clear robust
load(fullfile(params.sessionDir,b{runNum},'mc/motion_params.txt'));
robust = motion_params;
robust(:,1:3) = robust(:,1:3)*50;
%% Compare mri_robust_register and mcflirt
close all;

fullFigure;

subplot(2,2,1);
plot(robust);
ylim([-.5 .5]);
title('robust register','FontSize',20);
xlabel('TR','FontSize',20);
ylabel('movement (mm)','FontSize',20);

[h1,h2] = legend({'pitch' 'yaw' 'roll' 'x' 'y' 'z'},'FontSize',20);
set(h1,'Position',[0.49 0.5 .05 .1]);
hL=findobj(h2,'type','line');  % get the lines
set(hL,'linewidth',5);
hF=findobj(h2,'type','text');  % get the text
set(hF,'fontsize',20);

subplot(2,2,2);
plot(mcflirtOut);
ylim([-.5 .5]);
title('mcflirt','FontSize',20);
xlabel('TR','FontSize',20);
ylabel('movement (mm)','FontSize',20);


subplot(2,2,3)
x = detrend(robust,'constant');
pwelch(x,[],[],[],Fs);
xlim([0 700]);
ylim([-80 20]);
title('robust register','FontSize',20);
xlabel('Frequency (mHz)','FontSize',20);
ylabel('Power/Frequency (dB/Hz)','FontSize',20);

% x = robust;
% L = length(x); % Length of signal
% xdft = fft(x);
% Pxx = 1/(L*Fs)*abs(xdft(1:length(x)/2+1)).^2;
% freq = 0:Fs/L:Fs/2;
% plot(freq,10*log10(Pxx));
% xlabel('Hz'); ylabel('dB/Hz');

subplot(2,2,4);
x = detrend(mcflirtOut,'constant');
pwelch(x,[],[],[],Fs);
xlim([0 700]);
ylim([-80 20]);
title('mcflirt','FontSize',20);
xlabel('Frequency (mHz)','FontSize',20);
ylabel('Power/Frequency (dB/Hz)','FontSize',20);

cd(figDir);
savefigs('pdf');
close all;
% subplot(2,2,3);
% Fs = 199; % Sampling frequency
% x = detrend(pulse.pulse.data,'constant');
% %x = pulse.data;
% pwelch(x,[],[],[],Fs);
% xlim([0 2]);
% ylim([-20 100]);
% title('pulse');
%
% subplot(2,2,4);
% Fs = 199; % Sampling frequency
% x = detrend(pulse.pulse.data,'constant');
% %x = pulse.data;
% pwelch(x,[],[],[],Fs);
% xlim([0 2]);
% ylim([-20 100]);
% title('pulse');
%% Pulse data
fullFigure;
outDir = fullfile(params.sessionDir,b{runNum});
dicomDir = fullfile(params.sessionDir,'DICOMS',b{runNum});
pulsDir = fullfile(params.sessionDir,'PulseOx');
[pulsMatch] = matchPulsFile(dicomDir,pulsDir);
if ~isempty(pulsMatch)
    pulse = PulseResp(dicomDir,pulsMatch,outDir);
else
    pulse.all = [];
end
plot(pulse.pulse.data);
%% Plot Pulse data
close all;
fullFigure;
subplot(2,1,1)
x = detrend(robust,'constant');
pwelch(x,[],[],[],Fs);
xlim([0 700]);
ylim([-80 20]);
title('robust register','FontSize',20);
xlabel('Frequency (mHz)','FontSize',20);
ylabel('Power/Frequency (dB/Hz)','FontSize',20);

subplot(2,1,2)
x = detrend(pulse.pulse.data,'constant');
pwelch(x,[],[],[],pulse.pulse.sampR);
xlim([0 0.7]);
%ylim([0 100]);
title('Pulse Ox','FontSize',20);
xlabel('Frequency (Hz)','FontSize',20);
ylabel('Power/Frequency (dB/Hz)','FontSize',20);
cd(figDir);
savefigs('pdf');
close all;
%% mcflirt motion movie
yLimVals                    = [-1 1];
outDir                      = fullfile(params.sessionDir,b{runNum});
outFile                     = fullfile(outDir,'mcflirt');
outObj                      = VideoWriter(outFile);
outObj.FrameRate            = 30;
outObj.Quality              = 100;
open(outObj);
close all;
ih = figure;
for i = 1:size(mcflirtOut,1)
    hold off;
    plot(mcflirtOut);
    title('mcflirt','FontSize',20);
    xlim([0 size(mcflirtOut,1)]);
    xlabel('TR','FontSize',20);
    ylabel('movement (mm)','FontSize',20);
    axis square;
    hold on;
    plot([i-1 i-1],yLimVals,'r','LineWidth',2);
    ylim(yLimVals);
    frame                   = getframe(ih);
    writeVideo(outObj,frame);
end
close(outObj);
close(ih);
%% robust motion movie
yLimVals                    = [-1 1];
outDir                      = fullfile(params.sessionDir,b{runNum});
outFile                     = fullfile(outDir,'robust');
outObj                      = VideoWriter(outFile);
outObj.FrameRate            = 30;
outObj.Quality              = 100;
open(outObj);
close all;
ih = figure;
for i = 1:size(robust,1)
    hold off;
    plot(robust);
    title('robust','FontSize',20);
    xlim([0 size(robust,1)]);
    xlabel('TR','FontSize',20);
    ylabel('movement (mm)','FontSize',20);
    axis square;
    hold on;
    plot([i-1 i-1],yLimVals,'r','LineWidth',2);
    ylim(yLimVals);
    frame                   = getframe(ih);
    writeVideo(outObj,frame);
end
close(outObj);
close(ih);
%% mcflirt movie - brain
outDir                      = fullfile(params.sessionDir,b{runNum});
inFile                      = 'mcflirt';
outFile                     = fullfile(outDir,inFile);
outObj                      = VideoWriter(fullfile(outDir,'mcflirtBrain'));
outObj.FrameRate            = 30;
outObj.Quality              = 100;
% Convert to RAS
system(['mri_convert ' outFile '.nii.gz ' outFile '.RAS.nii.gz']);
% Load file
foo = load_nifti([outFile '.RAS.nii.gz']);
% Open video file
open(outObj);
% make video
ih = figure;
for i = 1:size(foo.vol,4);
    tmp                     = squeeze(foo.vol(size(foo.vol,1)/2,:,:,i));
    imagesc(squeeze(flipud(tmp')));
    colormap('gray');
    frame                   = getframe(ih);
    writeVideo(outObj,frame);
end
close(outObj);
close(ih);
%% robust movie - brain
outDir                      = fullfile(params.sessionDir,b{runNum});
inFile                      = 'rf';
outFile                     = fullfile(outDir,inFile);
outObj                      = VideoWriter(fullfile(outDir,'robustBrain'));
outObj.FrameRate            = 30;
outObj.Quality              = 100;
% Convert to RAS
system(['mri_convert ' outFile '.nii.gz ' outFile '.RAS.nii.gz']);
% Load file
foo = load_nifti([outFile '.RAS.nii.gz']);
% Open video file
open(outObj);
% make video
ih = figure;
for i = 1:size(foo.vol,4);
    tmp                     = squeeze(foo.vol(size(foo.vol,1)/2,:,:,i));
    imagesc(squeeze(flipud(tmp')));
    colormap('gray');
    frame                   = getframe(ih);
    writeVideo(outObj,frame);
end
close(outObj);
close(ih);
%% robust movie - brain
outDir                      = fullfile(params.sessionDir,b{runNum});
inFile                      = 'raw_f';
outFile                     = fullfile(outDir,inFile);
outObj                      = VideoWriter(fullfile(outDir,'rawBrain'));
outObj.FrameRate            = 30;
outObj.Quality              = 100;
% Convert to RAS
system(['mri_convert ' outFile '.nii.gz ' outFile '.RAS.nii.gz']);
% Load file
foo = load_nifti([outFile '.RAS.nii.gz']);
% Open video file
open(outObj);
% make video
ih = figure;
for i = 1:size(foo.vol,4);
    tmp                     = squeeze(foo.vol(size(foo.vol,1)/2,:,:,i));
    imagesc(squeeze(flipud(tmp')));
    colormap('gray');
    frame                   = getframe(ih);
    writeVideo(outObj,frame);
end
close(outObj);
close(ih);
%% Simulations
tmp                         = load_nifti(fullfile(params.sessionDir,'zeros.nii.gz'));
dims                        = size(tmp.vol);
x                           = 1:dims(1);
y                           = 1:dims(2);
z                           = 1:dims(3);
[xGrid,yGrid,zGrid]         = meshgrid(x,y,z);
dists                       = sqrt( (xGrid - dims(1)/2).^2 + (yGrid - dims(2)/2).^2 + ...
    (zGrid - dims(3)/2).^2 );
% Make a sphere
outFlat                     = (10^-6)*ones(dims(1)*dims(2)*dims(3),dims(4));
for i = 1:dims(4)
   outFlat(dists<20,i)        = 1;
end
out                         = tmp;
out.vol                     = reshape(outFlat,dims(1),dims(2),dims(3),dims(4));
save_nifti(out,fullfile(params.sessionDir,'sphere.nii.gz'));
% Make a hollow sphere
outFlat                     = zeros(dims(1)*dims(2)*dims(3),dims(4));
for i = 1:dims(4)
   outFlat(dists<20,i)        = 1000;
   outFlat(dists<10,i)        = 0;
end
out                         = tmp;
out.vol                     = reshape(outFlat,dims(1),dims(2),dims(3),dims(4));
save_nifti(out,fullfile(params.sessionDir,'hollow.nii.gz'));
%% Get motion params - robust register
params.despike                  = 0; % params.despike data
params.slicetiming              = 0; % do slice timing correction
params.refvol                   = 1; % reference volume = 1st TR
params.regFirst                 = 0; % register to the first run
runNum                          = 1;
b                               = find_bold(params.sessionDir);
% Correct motion using mri_robust_register
motion_slice_correction(params,runNum);
%% Get motion params - mcflirt
clear mcflirtOut
thisDir = fullfile(params.sessionDir);
mats = listdir(fullfile(thisDir,'sphere.10.mcflirt.nii.gz.mat++++/MAT*'),'files');
cd(fullfile(thisDir,'sphere.10.mcflirt.nii.gz.mat++++'));
mcflirtOut = nan(length(mats),6);
for i = 1:length(mats)
    [x,y,z,pitch,yaw,roll] = convertMAT2tranrot(mats{i});
    mcflirtOut(i,1) = pitch*50;
    mcflirtOut(i,2) = yaw*50;
    mcflirtOut(i,3) = roll*50;
    mcflirtOut(i,4) = x;
    mcflirtOut(i,5) = y;
    mcflirtOut(i,6) = z;
end
%%
close all;

fullFigure;

subplot(1,2,1);
plot(robust);
ylim([-.1 .1]);
title('robust register','FontSize',20);
xlabel('TR','FontSize',20);
ylabel('movement (mm)','FontSize',20);

[h1,h2] = legend({'pitch' 'yaw' 'roll' 'x' 'y' 'z'},'FontSize',20);
set(h1,'Position',[0.49 0.5 .05 .1]);
hL=findobj(h2,'type','line');  % get the lines
set(hL,'linewidth',5);
hF=findobj(h2,'type','text');  % get the text
set(hF,'fontsize',20);

subplot(1,2,2);
plot(mcflirtOut);
ylim([-.1 .1]);
title('mcflirt','FontSize',20);
xlabel('TR','FontSize',20);
ylabel('movement (mm)','FontSize',20);
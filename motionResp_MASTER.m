%% Script to investigate the influence of respiration on motion correction
%
%   Written by Andrew S Bock Oct 2016

%% set defaults
figDir                          = '/Users/abock/MOTION_CORRECTION_TECHDEV/figures';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/heteroPhantom/042016/';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/heteroPhantom/051016';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/heteroPhantom/052016';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/homoPhantom/042016';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/HCLV1001/8792BT';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/FBIRN/052516';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/FBIRN/100416';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3001/081916a';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3002/082616a';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3003/090216';
%params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3004/091916';
params.sessionDir               = '/Users/abock/MOTION_CORRECTION_TECHDEV/TOME_3005/092316';
params.despike                  = 0; % params.despike data
params.slicetiming              = 0; % do slice timing correction
params.refvol                   = 1; % reference volume = 1st TR
params.regFirst                 = 0; % register to the first run
runNum                          = 2;
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
ylim([-5 5]);
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
ylim([-5 5]);
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
%%
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
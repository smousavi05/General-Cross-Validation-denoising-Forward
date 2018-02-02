%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The following is a demo to run the modified version of GCV Denoising algorithm presented in:
%
%   Mousavi, S. M., and C. A. Langston (2017). Automatic Noise-Removal/
%   Signal-Removal Based on the General-Cross-Validation Thresholding in 
%   Synchrosqueezed domains, and its application on earthquake data, 
%   Geophysics.82(4), V211-V227 doi: 10.1190/geo2016-0433.1
%
%   Mousavi, S. M., and C. A. Langston, (2016a). Fast and novel microseismic 
%   detection using time-frequency analysis. SEG Technical Program Expanded 
%   Abstracts 2016: pp. 2632-2636. doi: 10.1190/segam2016-13262278.1                                                                                                                                                               
%
%  The modified version runs on contineous wavelet transform and just include
%  pre-processing and the main thresholding steps based on GCV and softThresholding.
%  Thus it is faster than the original algorithm. This algorithm does not rely on any
%  arrival time estimation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Parameters for Calculate the wavelet transform -
opt.type = 'morlet';         % Type od the mother wavelet used for CWT calculation: 'morlet' 'shannon' 'mhat' 'hhat'
opt.padtype = 'symmetric';   % padded via symmetrization
opt.rpadded = 1;
opt.nv = 32;                 % Number of voices. You can sellect 8,16,32,64.

% Guassian correction factor. This corrects the uncertainties for the 
% estimation of the guassianity and amplifies the pre-processing effects.
% It should be selected highh enough to remove the strong noises outside
% of the signal's frequency band but not too high to remove signal's energy. 
% value of 1.0 means no correction. 
opt.gc=3.5;

% processing long duration data is done in moving window fasion
 opt.wsiz = 120; % wondow size (sec), needs to be longer than the length of typical events that you have in your data
 
load('syntNoisy3_z.mat');
data.noisy = syntNoisy3_z;
data.t = linspace(0,(100),length(data.noisy));
data.dt = 0.1;

tic
dn =gcvThreshF(data,opt);
toc

figure, 
subplot(2,1,1),
plot(data.noisy, 'b');
grid on
title('Noisy Record ','Rotation',0,'FontSize',14);
xlabel({'Sample'},'FontSize',12); 
ylabel('Amplitude (count)','FontSize',12)


subplot(2,1,2),
plot(dn.xgcv, 'b');
grid on
title('Denoised Record ','Rotation',0,'FontSize',14);
xlabel({'Sample'},'FontSize',12); 
ylabel('Amplitude (count)','FontSize',12)


colormap('jet')


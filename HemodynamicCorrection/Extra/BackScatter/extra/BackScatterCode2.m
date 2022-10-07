%path and some parameters
addpath(genpath('Z:/GRABS_Data/Code/meso_processing-master'));
addpath(genpath('F:\Sweyta\GRABS_Data\BackScatter'));
rmpath(genpath('Z:\GRABS_Data\Code\meso_processing-master\chronux_2_12')); 
inputPath= 'F:\Sweyta\GRABS_Data\BackScatter';
outputPath='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\Figures\Hemodynamics\BackScatterOutput;
fsspike2=5000;
fsimaging=10;
pupilSR=10;
deterend.filtLen = 1000;
% visusal stim
visStimAn=1; % 1 if vis stim are presented, 0 if not
airpuffAn=0; % if airpuffs or electrical stim are presented, indicate 1;
visStimDur=2;
visStimITI=5;
%run params
minRunDuration=5; %in seconds
minSitDuration=10; %in seconds
ITITime=5; %  dead time after events in seconds
%analysis window for events (airpuff and visstim)
preEventWin=2;
postEventWin=5;
baselineWin=2; %use this period to calculate baseline during within trial normalization
% spike2 channels for bcmm
channels.BLUE = 1;%blue led
channels.UV = 2;%uv led
channels.GREEN=9;%green LED
channels.BS1=6; 
channels.BS2=10
channels.FRAMETICKS = 8;%green/uv mesocam triggers
channels.PHOTO_DIODE = 4;%visual stim trigger
channels.WHEEL = 5;
channels.PUPIL_CAMERA = 7;
channels.RED_MESO = 3;%red mesocam trigger


%% Method 1: Uncorrected load uncorrected dF/F

load(fullfile(outputPath,'dFOF_Uncorrected.mat'));
%get dFoF parcells
names=fieldnames(dFoF);
load('parcells_updated121519.mat'); parcells=parcells_new;
for i = 1:length(names)
    disp(['Parcellating' names{i}]);
    dFoF_parcells.(names{i}) = pixels2rois(dFoF.(names{i}), parcells);
end

%% Method 2: Hadas SVD against UV

load(fullfile(outputPath,'dFOF_SVD_UV_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_SVD_UV,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_SVD_UV_parcells.(names{i}) = pixels2rois(dFoF_SVD_UV.(names{i}), parcells);
end

%% Method 3: Hadas SVD against BS1

load(fullfile(outputPath,'dFOF_SVD_BS1_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_SVD_BS1,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_SVD_BS1_parcells.(names{i}) = pixels2rois(dFoF_SVD_BS1.(names{i}), parcells);
end

%%Method 4: pixelwise regression against UV
load(fullfile(outputPath,'dFOF_PixUV_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_PixUV,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_PixUV_parcells.(names{i}) = pixels2rois(dFoF_PixUV.(names{i}), parcells);
end

%% Method 5: pixelwise regression against BS1

load(fullfile(outputPath,'dFOF_PixBS1_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_PixBS1,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_PixBS1_parcells.(names{i}) = pixels2rois(dFoF_PixBS1.(names{i}), parcells);
end


%% Method 6: pixelwise regression against BS1 and BS2
load(fullfile(outputPath,'dFOF_PixBS12_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_PixBS12,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_PixBS12_parcells.(names{i}) = pixels2rois(dFoF_PixBS12.(names{i}), parcells);
end

%% Method 7: Beer Lambert model (Waters and Hillman paper)
load(fullfile(outputPath,'dFOF_BLam_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_BLam,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_BLam_parcells.(names{i}) = pixels2rois(dFoF_BLam.(names{i}), parcells);
end

%% Method 8: Hadas SVD model against BS2
load(fullfile(outputPath,'dFOF_SVD_BS2_nosmooth.mat'));
names=fieldnames(dFoF_SVD_BS2);
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_BLam,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_SVD_BS2_parcells.(names{i}) = pixels2rois(dFoF_SVD_BS2.(names{i}), parcells);
end

%% Method 9: Pixelwise regress against BS2
%load(fullfile(outputPath,'dFOF_PixBS2_nosmooth.mat'));
names=fieldnames(dFoF_PixBS2);
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_BLam,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_PixBS2_parcells.(names{i}) = pixels2rois(dFoF_PixBS2.(names{i}), parcells);
end

%% calculate R2 comparing traces for V1, S1, and M2 parcells across methods
%%comapre all to the SVD against UV method 
Idx.visual=2:2:16; %visual parcells
Idx.PosteriorParietal=32;% PPC
Idx.RSL_AgL=[18,20];%retrosplenial cortex
Idx.Somat=34:2:48;%somatosensory areas
Idx.FrontalMotor=50:2:52;%frontal-moto area
Idx.auditory=[28,30];%auditory areas
CombinedParcellIdx=[Idx.visual,Idx.RSL_AgL,Idx.auditory,Idx.PosteriorParietal,Idx.Somat,Idx.FrontalMotor];

names={'blue','green'};
for i=1:length(names)
    for j=1:length(parcells.names)
        R2_par.(names{i})(j,1)=corr(dFoF_parcells.(names{i})(j,:)',dFoF_SVD_UV_parcells.(names{i})(j,:)').^2; %method 1 uncorrected
        R2_par.(names{i})(j,2)=corr(dFoF_SVD_UV_parcells.(names{i})(j,:)',dFoF_SVD_UV_parcells.(names{i})(j,:)').^2; %method 2 SVD UV
        R2_par.(names{i})(j,3)=corr(dFoF_SVD_BS1_parcells.(names{i})(j,:)',dFoF_SVD_UV_parcells.(names{i})(j,:)').^2; %method 3 SVD BS1
        R2_par.(names{i})(j,4)=corr(dFoF_PixUV_parcells.(names{i})(j,:)',dFoF_SVD_UV_parcells.(names{i})(j,:)').^2; %method 4 pixelwise uv regression
        R2_par.(names{i})(j,5)=corr(dFoF_PixBS1_parcells.(names{i})(j,:)',dFoF_SVD_UV_parcells.(names{i})(j,:)').^2; %method 5 pixelwise BS1 regression
        R2_par.(names{i})(j,6)=corr(dFoF_PixBS12_parcells.(names{i})(j,:)',dFoF_SVD_UV_parcells.(names{i})(j,:)').^2; %method 6 pixelwise BS1 and BS2 regression
        R2_par.(names{i})(j,7)=corr(dFoF_BLam_parcells.(names{i})(j,:)',dFoF_SVD_UV_parcells.(names{i})(j,:)').^2; %method 7 beer lambert model 
    end    
    
end
figure; suptitle('R2-comparedwithSVDUV-GRABs');
subplot(7,1,1);bar(R2_par.blue(CombinedParcellIdx,1)); ylim([0 1]); title('NonCorrected');
subplot(7,1,2);bar(R2_par.blue(CombinedParcellIdx,2)); ylim([0 1]); title('SVD UV')
subplot(7,1,3);bar(R2_par.blue(CombinedParcellIdx,3)); ylim([0 1]); title('SVD BS1')
subplot(7,1,4);bar(R2_par.blue(CombinedParcellIdx,4)); ylim([0 1]); title('Pixelwise Regress UV')
subplot(7,1,5);bar(R2_par.blue(CombinedParcellIdx,5)); ylim([0 1]); title('Pixelwise Regress BS1')
subplot(7,1,6);bar(R2_par.blue(CombinedParcellIdx,6)); ylim([0 1]); title('Pixelwise Regress BS1 and 2')
subplot(7,1,7);bar(R2_par.blue(CombinedParcellIdx,7)); ylim([0 1]); title('Beer Lamber-Simplified');xticks(1:length(CombinedParcellIdx));xticklabels(parcells.names(CombinedParcellIdx));xtickangle(90);box off; set(gca,'TickDir','out')

figure; suptitle('R2-comparedwithSVDUV-RCaMP');
subplot(7,1,1);bar(R2_par.green(CombinedParcellIdx,1)); ylim([0 1]); title('NonCorrected');
subplot(7,1,2);bar(R2_par.green(CombinedParcellIdx,2)); ylim([0 1]); title('SVD UV')
subplot(7,1,3);bar(R2_par.green(CombinedParcellIdx,3)); ylim([0 1]); title('SVD BS1')
subplot(7,1,4);bar(R2_par.green(CombinedParcellIdx,4)); ylim([0 1]); title('Pixelwise Regress UV')
subplot(7,1,5);bar(R2_par.green(CombinedParcellIdx,5)); ylim([0 1]); title('Pixelwise Regress BS1')
subplot(7,1,6);bar(R2_par.green(CombinedParcellIdx,6)); ylim([0 1]); title('Pixelwise Regress BS1 and 2')
subplot(7,1,7);bar(R2_par.green(CombinedParcellIdx,7)); ylim([0 1]); title('Beer Lamber-Simplified');xticks(1:length(CombinedParcellIdx));xticklabels(parcells.names(CombinedParcellIdx));xtickangle(90);box off; set(gca,'TickDir','out')

figure; suptitle('TraceV1-GRABs');
ax1=subplot(7,1,1);plot(dFoF_parcells.blue(2,1000:3000)); title('NonCorrected');box off; ylim([-0.005 0.005]); set(gca,'TickDir','out')
ax2=subplot(7,1,2);plot(dFoF_SVD_UV_parcells.blue(2,1000:3000)); title('SVD UV');box off; ylim([-0.005 0.005]);set(gca,'TickDir','out')
ax3=subplot(7,1,3);plot(dFoF_SVD_BS1_parcells.blue(2,1000:3000)); title('SVD BS1');box off; ylim([-0.005 0.005]);set(gca,'TickDir','out')
ax4=subplot(7,1,4);plot(dFoF_PixUV_parcells.blue(2,1000:3000));title('Pixelwise Regress UV');ylim([-0.005 0.005]);box off; set(gca,'TickDir','out')
ax5=subplot(7,1,5);plot(dFoF_PixBS1_parcells.blue(2,1000:3000)); title('Pixelwise Regress BS1');ylim([-0.005 0.005]);box off; set(gca,'TickDir','out')
ax6=subplot(7,1,6);plot(dFoF_PixBS12_parcells.blue(2,1000:3000));title('Pixelwise Regress BS1 and 2');ylim([-0.005 0.005]);box off; set(gca,'TickDir','out')
ax7=subplot(7,1,7);plot(dFoF_BLam_parcells.blue(2,1000:3000)); title('Beer Lamber-Simplified');box off; ylim([-0.005 0.005]);set(gca,'TickDir','out')
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'xy')
figure; suptitle('TraceV1-RCaMP');
ax1=subplot(7,1,1);plot(dFoF_parcells.green(2,1000:3000)); title('NonCorrected');box off; ylim([-0.008 0.008]);set(gca,'TickDir','out')
ax2=subplot(7,1,2);plot(dFoF_SVD_UV_parcells.green(2,1000:3000)); title('SVD UV');box off; ylim([-0.008 0.008]);set(gca,'TickDir','out')
ax3=subplot(7,1,3);plot(dFoF_SVD_BS1_parcells.green(2,1000:3000)); title('SVD BS1');box off; ylim([-0.008 0.008]);set(gca,'TickDir','out')
ax4=subplot(7,1,4);plot(dFoF_PixUV_parcells.green(2,1000:3000));title('Pixelwise Regress UV');ylim([-0.008 0.008]);box off; set(gca,'TickDir','out')
ax5=subplot(7,1,5);plot(dFoF_PixBS1_parcells.green(2,1000:3000)); title('Pixelwise Regress BS1');ylim([-0.008 0.008]);box off; set(gca,'TickDir','out')
ax6=subplot(7,1,6);plot(dFoF_PixBS12_parcells.green(2,1000:3000));title('Pixelwise Regress BS1 and 2');ylim([-0.008 0.008]);box off; set(gca,'TickDir','out')
ax7=subplot(7,1,7);plot(dFoF_BLam_parcells.green(2,1000:3000)); title('Beer Lamber-Simplified');box off; ylim([-0.008 0.008]);set(gca,'TickDir','out')
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'xy')
figure; suptitle('TraceM2-GRABs');
ax1=subplot(7,1,1);plot(dFoF_parcells.blue(52,1000:3000)); title('NonCorrected');box off; ylim([-0.01 0.01]); set(gca,'TickDir','out')
ax2=subplot(7,1,2);plot(dFoF_SVD_UV_parcells.blue(52,1000:3000)); title('SVD UV');box off; ylim([-0.01 0.01]);set(gca,'TickDir','out')
ax3=subplot(7,1,3);plot(dFoF_SVD_BS1_parcells.blue(52,1000:3000)); title('SVD BS1');box off; ylim([-0.01 0.01]);set(gca,'TickDir','out')
ax4=subplot(7,1,4);plot(dFoF_PixUV_parcells.blue(52,1000:3000));title('Pixelwise Regress UV');ylim([-0.01 0.01]);box off; set(gca,'TickDir','out')
ax5=subplot(7,1,5);plot(dFoF_PixBS1_parcells.blue(52,1000:3000)); title('Pixelwise Regress BS1');ylim([-0.01 0.01]);box off; set(gca,'TickDir','out')
ax6=subplot(7,1,6);plot(dFoF_PixBS12_parcells.blue(52,1000:3000));title('Pixelwise Regress BS1 and 2');ylim([-0.01 0.01]);box off; set(gca,'TickDir','out')
ax7=subplot(7,1,7);plot(dFoF_BLam_parcells.blue(52,1000:3000)); title('Beer Lamber-Simplified');box off; ylim([-0.01 0.01]);set(gca,'TickDir','out')
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'xy')
figure; suptitle('TraceM2-RCaMP');
ax1=subplot(7,1,1);plot(dFoF_parcells.green(52,1000:3000)); title('NonCorrected');box off; ylim([-0.01 0.01]);set(gca,'TickDir','out')
ax2=subplot(7,1,2);plot(dFoF_SVD_UV_parcells.green(52,1000:3000)); title('SVD UV');box off; ylim([-0.01 0.01]);set(gca,'TickDir','out')
ax3=subplot(7,1,3);plot(dFoF_SVD_BS1_parcells.green(52,1000:3000)); title('SVD BS1');box off; ylim([-0.01 0.01]);set(gca,'TickDir','out')
ax4=subplot(7,1,4);plot(dFoF_PixUV_parcells.green(52,1000:3000));title('Pixelwise Regress UV');ylim([-0.01 0.01]);box off; set(gca,'TickDir','out')
ax5=subplot(7,1,5);plot(dFoF_PixBS1_parcells.green(52,1000:3000)); title('Pixelwise Regress BS1');ylim([-0.01 0.01]);box off; set(gca,'TickDir','out')
ax6=subplot(7,1,6);plot(dFoF_PixBS12_parcells.green(52,1000:3000));title('Pixelwise Regress BS1 and 2');ylim([-0.01 0.01]);box off; set(gca,'TickDir','out')
ax7=subplot(7,1,7);plot(dFoF_BLam_parcells.green(52,1000:3000)); title('Beer Lamber-Simplified');box off; ylim([-0.01 0.01]);set(gca,'TickDir','out')
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'xy')


%%calculate R2 comparing traces for V1, S1, and M2 parcells across method for the whole brain

%%comapre all to the SVD against UV method 
load('brainMask.mat','mask');
maskinds = find(mask);
R=256; C=256; 
names={'blue','green'};
for i=1:length(names)
    R2.(names{i})=zeros(R*C,7);
    for j=1:length(maskinds)
        idx=maskinds(j); 
        R2.(names{i})(idx,1)=corr(dFoF.(names{i})(idx,:)',dFoF_SVD_UV.(names{i})(idx,:)').^2; %method 1 uncorrected
        R2.(names{i})(idx,2)=corr(dFoF_SVD_UV.(names{i})(idx,:)',dFoF_SVD_UV.(names{i})(idx,:)').^2; %method 2 SVD UV
        R2.(names{i})(idx,3)=corr(dFoF_SVD_BS1.(names{i})(idx,:)',dFoF_SVD_UV.(names{i})(idx,:)').^2; %method 3 SVD BS1
        R2.(names{i})(idx,4)=corr(dFoF_PixUV.(names{i})(idx,:)',dFoF_SVD_UV.(names{i})(idx,:)').^2; %method 4 pixelwise uv regression
        R2.(names{i})(idx,5)=corr(dFoF_PixBS1.(names{i})(idx,:)',dFoF_SVD_UV.(names{i})(idx,:)').^2; %method 5 pixelwise BS1 regression
        R2.(names{i})(idx,6)=corr(dFoF_PixBS12.(names{i})(idx,:)',dFoF_SVD_UV.(names{i})(idx,:)').^2; %method 6 pixelwise BS1 and BS2 regression
        R2.(names{i})(idx,7)=corr(dFoF_BLam.(names{i})(idx,:)',dFoF_SVD_UV.(names{i})(idx,:)').^2; %method 7 beer lambert model 
    end    
    
end

figure; suptitle('R2-comparedwithSVDUV');
subplot(6,2,1);imagesc(reshape(R2.blue(:,1),R,C)); caxis([0 0.5]);title('GRABs-NonCorrected');
subplot(6,2,3);imagesc(reshape(R2.blue(:,3),R,C)); caxis([0 0.5]); title('SVD BS1')
subplot(6,2,5);imagesc(reshape(R2.blue(:,4),R,C)); caxis([0 0.5]); title('Pixelwise Regress UV')
subplot(6,2,7);imagesc(reshape(R2.blue(:,5),R,C)); caxis([0 0.5]); title('Pixelwise Regress BS1')
subplot(6,2,9);imagesc(reshape(R2.blue(:,6),R,C)); caxis([0 0.5]); title('Pixelwise Regress BS1 and 2')
subplot(6,2,11);imagesc(reshape(R2.blue(:,7),R,C)); caxis([0 0.5]); title('Beer Lamber-Simplified');

subplot(6,2,2);imagesc(reshape(R2.green(:,1),R,C)); caxis([0 0.5]); title('RCaMP=NonCorrected');
subplot(6,2,4);imagesc(reshape(R2.green(:,3),R,C)); caxis([0 0.5]); title('SVD BS1')
subplot(6,2,6);imagesc(reshape(R2.green(:,4),R,C)); caxis([0 0.5]); title('Pixelwise Regress UV')
subplot(6,2,8);imagesc(reshape(R2.green(:,5),R,C)); caxis([0 0.5]); title('Pixelwise Regress BS1')
subplot(6,2,10);imagesc(reshape(R2.green(:,6),R,C)); caxis([0 0.5]); title('Pixelwise Regress BS1 and 2')
subplot(6,2,12);imagesc(reshape(R2.green(:,7),R,C)); caxis([0 0.5]); title('Beer Lamber-Simplified');

%plot traces and get R2 for two pixels (one in visual and one in frontal) 
VisualPixel_2d=[58,188]; 
VisPix=sub2ind([R, C],VisualPixel_2d(1),VisualPixel_2d(2));

FrontalPixel_2d=[166,163]; 
FronPix=sub2ind([R, C],FrontalPixel_2d(1),FrontalPixel_2d(2));


%% plot trial averaged activity around visual stimulus in V1, S1, M2

%% load and process spike2 , get timestamps
smrxfilelist = (dir(fullfile(inputPath, '*.smrx')));
dataSmrxFile=fullfile(inputPath,smrxfilelist.name);
cedpath = '../pre_processing_scripts/utils/CEDS64ML';

display(strcat('loading spike 2 file: ',dataSmrxFile));
if exist(fullfile(outputPath, 'smrx_signals.mat'),'file')
    load(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx');
else
    [timing, channels_data] = process_spike2_two_cams_BS(cedpath, outputPath,dataSmrxFile, fsspike2,channels,minRunDuration,minSitDuration,ITITime,visStimDur,visStimITI,pupilSR);
    
    % get timestamps for all relevant events
    timestamps.timaging = (timing.bluestart+timing.blueend)/2;    
    timestamps.timaging=timestamps.timaging(deterend.filtLen/2:end);%remove the first filtLen/2 timestamps because there is a filtering artifact and those samples are removed in the detrending step
    timestamps.timaging=timestamps.timaging(1:size(dFoF.blue,2));%remove excess timestamps if you don't have corresponding frames
    timestamps.visStim=timing.stimstart;
    timestamps.wheelOn=timing.wheelOn(timing.wheelOn>(timestamps.timaging(1)+ minSitDuration) & timing.wheelOn<(timestamps.timaging(end)-minRunDuration));%remove events if they occur outside imaging period
    
    %% if you want to analyze visual stim
    contrasts=100;
    visStimidx= find(timing.stimstart>(timestamps.timaging(1)+preEventWin) & timing.stimstart<(timestamps.timaging(end)-postEventWin));
    timestamps.visStim =timing.stimstart(visStimidx);%remove events if they occur outside of imaging period
    contrasts=repmat(contrasts,1,length(visStimidx));
    uniqContrasts=unique(contrasts);
    for i=1:length(uniqContrasts)
        contrast_Idx{uniqContrasts(i)}=find(contrasts==uniqContrasts(i));
        timestamps.contrasts{uniqContrasts(i)}=timestamps.visStim(contrast_Idx{uniqContrasts(i)});
    end
    save(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx');
end

%% draw raw traces in visual parcell along with vis timestamps 
timeframes=1300:2000;
imaging_timestamps=timestamps.timaging(timeframes);
brainMapPlotTime=198.3; 
visual_times=timestamps.visStim(timestamps.visStim>imaging_timestamps(1) & timestamps.visStim<imaging_timestamps(end));
figure; suptitle('TraceV1-GRABs');
ax1=subplot(7,1,1);plot(imaging_timestamps,dFoF_parcells.blue(2,timeframes)); title('NonCorrected');box off; axis off; ylim([-0.005 0.005]);ylimits=ylim; set(gca,'TickDir','out'); hold on;
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
 plot([brainMapPlotTime,brainMapPlotTime], ylimits, 'm');
ax2=subplot(7,1,2);plot(imaging_timestamps,dFoF_SVD_UV_parcells.blue(2,timeframes)); title('SVD UV');box off; axis off; ylim([-0.005 0.005]);set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax3=subplot(7,1,3);plot(imaging_timestamps,dFoF_SVD_BS1_parcells.blue(2,timeframes)); title('SVD BS1');box off; axis off; ylim([-0.005 0.005]);set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax4=subplot(7,1,4);plot(imaging_timestamps,dFoF_PixUV_parcells.blue(2,timeframes));title('Pixelwise Regress UV');ylim([-0.005 0.005]);box off;axis off; set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax5=subplot(7,1,5);plot(imaging_timestamps,dFoF_PixBS1_parcells.blue(2,timeframes)); title('Pixelwise Regress BS1');ylim([-0.005 0.005]);box off; axis off; set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax6=subplot(7,1,6);plot(imaging_timestamps,dFoF_PixBS12_parcells.blue(2,timeframes));title('Pixelwise Regress BS1 and 2');ylim([-0.005 0.005]);box off;axis off;  set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax7=subplot(7,1,7);plot(imaging_timestamps,dFoF_BLam_parcells.blue(2,timeframes)); title('Beer Lamber-Simplified');box off; ylim([-0.005 0.005]);set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'xy')

figure; suptitle('TraceV1-RCaMP');
ax1=subplot(7,1,1);plot(imaging_timestamps,dFoF_parcells.green(2,timeframes)); title('NonCorrected');box off; axis off; ylim([-0.005 0.007]);ylimits=ylim; set(gca,'TickDir','out'); hold on;
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
 plot([brainMapPlotTime,brainMapPlotTime], ylimits, 'm');
ax2=subplot(7,1,2);plot(imaging_timestamps,dFoF_SVD_UV_parcells.green(2,timeframes)); title('SVD UV');box off; axis off; ylim([-0.005 0.007]);set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax3=subplot(7,1,3);plot(imaging_timestamps,dFoF_SVD_BS1_parcells.green(2,timeframes)); title('SVD BS1');box off; axis off; ylim([-0.005 0.007]);set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax4=subplot(7,1,4);plot(imaging_timestamps,dFoF_PixUV_parcells.green(2,timeframes));title('Pixelwise Regress UV');ylim([-0.005 0.007]);box off; axis off; set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax5=subplot(7,1,5);plot(imaging_timestamps,dFoF_PixBS1_parcells.green(2,timeframes)); title('Pixelwise Regress BS1');ylim([-0.005 0.007]);box off;axis off; set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax6=subplot(7,1,6);plot(imaging_timestamps,dFoF_PixBS12_parcells.green(2,timeframes));title('Pixelwise Regress BS1 and 2');ylim([-0.005 0.007]);box off;axis off;  set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
ax7=subplot(7,1,7);plot(imaging_timestamps,dFoF_BLam_parcells.green(2,timeframes)); title('Beer Lamber-Simplified');box off; ylim([-0.005 0.007]);set(gca,'TickDir','out'); hold on; 
for tt=1:length(visual_times)
    plot([visual_times(tt),visual_times(tt)], ylimits, 'k');
end
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'xy')

% also draw brain maps at one time point specified above 
f = @(x) imgaussfilt(x, 2);
[~,timepoint] =(min(abs(timestamps.timaging-brainMapPlotTime))); 
figure; suptitle('df/F');
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF.blue(maskinds,timepoint); 
subplot(7,2,1);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.005 0.005]); title('GRABs-Uncorrected');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_SVD_UV.blue(maskinds,timepoint); 
subplot(7,2,3);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.005 0.005]);title('SVD-UV');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_SVD_BS1.blue(maskinds,timepoint); 
subplot(7,2,5);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.005 0.005]); title('SVD-BS1');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_PixUV.blue(maskinds,timepoint); 
subplot(7,2,7);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.005 0.005]);title('Pixelwise Regress UV');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_PixBS1.blue(maskinds,timepoint); 
subplot(7,2,9);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.005 0.005]); title('Pixelwise Regress BS1');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_PixBS12.blue(maskinds,timepoint); 
subplot(7,2,11);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2)); caxis([-0.005 0.005]); title('Pixelwise Regress BS1 and 2');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_BLam.blue(maskinds,timepoint); 
subplot(7,2,13);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2)); caxis([-0.005 0.005]); title('Beer Lamber-Simplified');box off; axis off

toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF.green(maskinds,timepoint); 
subplot(7,2,2);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.01 0.01]); title('RCaMP-Uncorrected');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_SVD_UV.green(maskinds,timepoint); 
subplot(7,2,4);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.01 0.01]);title('SVD-UV');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_SVD_BS1.green(maskinds,timepoint); 
subplot(7,2,6);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.01 0.01]); title('SVD-BS1');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_PixUV.green(maskinds,timepoint); 
subplot(7,2,8);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.01 0.01]); title('Pixelwise Regress UV');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_PixBS1.green(maskinds,timepoint); 
subplot(7,2,10);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2));caxis([-0.01 0.01]);  title('Pixelwise Regress BS1');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_PixBS12.green(maskinds,timepoint); 
subplot(7,2,12);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2)); caxis([-0.01 0.01]); title('Pixelwise Regress BS1 and 2');box off; axis off
toPlotImage=nan(R*C,1); toPlotImage(maskinds)=dFoF_BLam.green(maskinds,timepoint); 
subplot(7,2,14);tmp=reshape(toPlotImage,R,C);I2 = imrotate(roifilt2(tmp,mask,f),180);imagesc(I2,'AlphaData',~isnan(I2)); caxis([-0.01 0.01]);  title('Beer Lamber-Simplified');box off; axis off



%% plot visual trial activity method 1
h1=figure; title('V1-GRABs');hold on; h2=figure; title('V1-RCaMP');hold on; 
h3=figure; title('M2-GRABs');hold on; h4=figure; title('M2-RCaMP');hold on; 
V1Idx=2; 
S1bIdx=34; 
M2Idx=52; 

V1data.blue=dFoF_parcells.blue(V1Idx,:);
S1data.blue=dFoF_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_parcells.blue(M2Idx,:);
V1data.green=dFoF_parcells.green(V1Idx,:);
S1data.green=dFoF_parcells.green(S1bIdx,:);
M2data.green=dFoF_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax1=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax2=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax3=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax4=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 


%% plot visual trial activity method 2
V1data.blue=dFoF_SVD_UV_parcells.blue(V1Idx,:);
S1data.blue=dFoF_SVD_UV_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_SVD_UV_parcells.blue(M2Idx,:);
V1data.green=dFoF_SVD_UV_parcells.green(V1Idx,:);
S1data.green=dFoF_SVD_UV_parcells.green(S1bIdx,:);
M2data.green=dFoF_SVD_UV_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax5=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax6=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax7=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax8=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 


%% plot visual trial activity method 3
V1data.blue=dFoF_SVD_BS1_parcells.blue(V1Idx,:);
S1data.blue=dFoF_SVD_BS1_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_SVD_BS1_parcells.blue(M2Idx,:);
V1data.green=dFoF_SVD_BS1_parcells.green(V1Idx,:);
S1data.green=dFoF_SVD_BS1_parcells.green(S1bIdx,:);
M2data.green=dFoF_SVD_BS1_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax9=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax10=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax11=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax12=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 

%% plot visual trial activity method 4
V1data.blue=dFoF_PixUV_parcells.blue(V1Idx,:);
S1data.blue=dFoF_PixUV_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_PixUV_parcells.blue(M2Idx,:);
V1data.green=dFoF_PixUV_parcells.green(V1Idx,:);
S1data.green=dFoF_PixUV_parcells.green(S1bIdx,:);
M2data.green=dFoF_PixUV_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax13=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax14=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax15=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax16=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 


%% plot visual trial activity method 5
V1data.blue=dFoF_PixBS1_parcells.blue(V1Idx,:);
S1data.blue=dFoF_PixBS1_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_PixBS1_parcells.blue(M2Idx,:);
V1data.green=dFoF_PixBS1_parcells.green(V1Idx,:);
S1data.green=dFoF_PixBS1_parcells.green(S1bIdx,:);
M2data.green=dFoF_PixBS1_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax17=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax18=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax19=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax20=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 

%% plot visual trial activity method 6
V1data.blue=dFoF_PixBS12_parcells.blue(V1Idx,:);
S1data.blue=dFoF_PixBS12_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_PixBS12_parcells.blue(M2Idx,:);
V1data.green=dFoF_PixBS12_parcells.green(V1Idx,:);
S1data.green=dFoF_PixBS12_parcells.green(S1bIdx,:);
M2data.green=dFoF_PixBS12_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax21=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax22=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax23=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax24=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 

%% plot visual trial activity method 7
V1data.blue=dFoF_BLam_parcells.blue(V1Idx,:);
S1data.blue=dFoF_BLam_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_BLam_parcells.blue(M2Idx,:);
V1data.green=dFoF_BLam_parcells.green(V1Idx,:);
S1data.green=dFoF_BLam_parcells.green(S1bIdx,:);
M2data.green=dFoF_BLam_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax25=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; legend([ax1,ax5,ax9,ax13,ax17,ax21,ax25],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h2); ax26=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
legend([ax2,ax6,ax10,ax14,ax18,ax22,ax26],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h3); ax27=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
legend([ax3,ax7,ax11,ax15,ax19,ax23,ax27],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h4); ax28=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 
legend([ax4,ax8,ax12,ax16,ax20,ax24,ax28],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 


%% plot visual trial activity method 8
V1data.blue=dFoF_SVD_BS2_parcells.blue(V1Idx,:);
S1data.blue=dFoF_SVD_BS2_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_SVD_BS2_parcells.blue(M2Idx,:);
V1data.green=dFoF_SVD_BS2_parcells.green(V1Idx,:);
S1data.green=dFoF_SVD_BS2_parcells.green(S1bIdx,:);
M2data.green=dFoF_SVD_BS2_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 
set(0,'CurrentFigure',h1); ax29=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off;hold on; 
set(0,'CurrentFigure',h2); ax30=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax31=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax32=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 

%% plot visual trial activity method 9
V1data.blue=dFoF_PixBS2_parcells.blue(V1Idx,:);
S1data.blue=dFoF_PixBS2_parcells.blue(S1bIdx,:);
M2data.blue=dFoF_PixBS2_parcells.blue(M2Idx,:);
V1data.green=dFoF_PixBS2_parcells.green(V1Idx,:);
S1data.green=dFoF_PixBS2_parcells.green(S1bIdx,:);
M2data.green=dFoF_PixBS2_parcells.green(M2Idx,:);

%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.blue,timeStamps]= TrialTimeArrangeDff(S1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.blue,timeStamps]= TrialTimeArrangeDff(M2data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[S1trial.green,timeStamps]= TrialTimeArrangeDff(S1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[M2trial.green,timeStamps]= TrialTimeArrangeDff(M2data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%Zscore and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.blue(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.blue(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.blue-preEventMean)./preEventStd; 
M2trial_aveZ.blue=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZ.green=nanmean(Zscoredtrial,2); 
preEventMean=nanmean(M2trial.green(1:preEventWin*fsimaging,:)); 
preEventStd=nanstd(M2trial.green(1:preEventWin*fsimaging,:)); 
Zscoredtrial=(M2trial.green-preEventMean)./preEventStd; 
M2trial_aveZ.green=nanmean(Zscoredtrial,2); 

set(0,'CurrentFigure',h1); ax33=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.blue);box off; hold on; 
legend([ax1,ax5,ax9,ax13,ax17,ax21,ax25,ax29,ax33],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified','SpatialSVD-BS2','PixelwiseRegress-BS2'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h2); ax34=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZ.green);box off; hold on; 
legend([ax2,ax6,ax10,ax14,ax18,ax22,ax26,ax30,ax34],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified','SpatialSVD-BS2','PixelwiseRegress-BS2'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h3); ax35=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.blue);box off; hold on; 
legend([ax3,ax7,ax11,ax15,ax19,ax23,ax27,ax31,ax35],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified','SpatialSVD-BS2','PixelwiseRegress-BS2'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h4); ax36=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',M2trial_aveZ.green);box off; hold on; 
legend([ax4,ax8,ax12,ax16,ax20,ax24,ax28,ax32,ax36],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1','PixelwiseRegress-BS1-BS2','Beer-LamberSimplified','SpatialSVD-BS2','PixelwiseRegress-BS2'}); legend boxoff; 
ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 

saveas(h1,fullfile(outputPath,'VisualStimAverage-V1-GRABs')); saveas(h1,fullfile(outputPath,'VisualStimAverage-V1-GRABs.pdf')); 
saveas(h2,fullfile(outputPath,'VisualStimAverage-V1-RCaMP')); saveas(h2,fullfile(outputPath,'VisualStimAverage-V1-RCaMP.pdf')); 
saveas(h3,fullfile(outputPath,'VisualStimAverage-M2-GRABs')); saveas(h3,fullfile(outputPath,'VisualStimAverage-M2-GRABs.pdf')); 
saveas(h4,fullfile(outputPath,'VisualStimAverage-M2-RCaMP')); saveas(h4,fullfile(outputPath,'VisualStimAverage-M2=RCaMP.pdf')); 


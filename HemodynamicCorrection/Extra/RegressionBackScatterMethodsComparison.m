%path and some parameters
addpath(genpath('F:\Sweyta\GRABS_Data\FinalPaperCodes\PreprocessingCode'));
addpath(genpath('F:\Sweyta\GRABS_Data\BackScatter'));
inputPath= 'Z:\GRABS_Data\Analyzed_SVDMethodPatch14\Backscatter';
outputPath='F:\GRABS_Data\Analyzed_SVDMethodPatch14\Figures\Revisions2\Hemodynamics\RegressionBackscatterMethodComparison';
if ~exist(outputPath,'dir'),mkdir(outputPath); end 
%% parameters 
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
channels.BS2=10;
channels.FRAMETICKS = 8;%green/uv mesocam triggers
channels.PHOTO_DIODE = 4;%visual stim trigger
channels.WHEEL = 5;
channels.PUPIL_CAMERA = 7;
channels.RED_MESO = 3;%red mesocam trigger

%% Method 1: Uncorrected load uncorrected dF/F
load(fullfile(inputPath,'dFOF_Uncorrected.mat'));
%get dFoF parcells
names=fieldnames(dFoF);
load('parcells_updated121519.mat'); parcells=parcells_new;
for i = 1:length(names)
    disp(['Parcellating' names{i}]);
    dFoF_parcells.(names{i}) = pixels2rois(dFoF.(names{i}), parcells);
end

%% Method 2: Hadas SVD against UV
load(fullfile(inputPath,'dFOF_SVD_UV_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_SVD_UV,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_SVD_UV_parcells.(names{i}) = pixels2rois(dFoF_SVD_UV.(names{i}), parcells);
end

%% Method 3: Hadas SVD against BS1
load(fullfile(inputPath,'dFOF_SVD_BS1_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_SVD_BS1,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_SVD_BS1_parcells.(names{i}) = pixels2rois(dFoF_SVD_BS1.(names{i}), parcells);
end

%% %Method 4: pixelwise regression against UV
load(fullfile(inputPath,'dFOF_PixUV_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_PixUV,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_PixUV_parcells.(names{i}) = pixels2rois(dFoF_PixUV.(names{i}), parcells);
end

%% Method 5: pixelwise regression against BS1
load(fullfile(inputPath,'dFOF_PixBS1_nosmooth.mat'));
%get dFoF parcells
for i = 1:length(names)
    if ~isfield(dFoF_PixBS1,names{i}),continue,end 
    disp(['Parcellating' names{i}]);
    dFoF_PixBS1_parcells.(names{i}) = pixels2rois(dFoF_PixBS1.(names{i}), parcells);
end

%% load and process spike2,get timestamps
smrxfilelist = (dir(fullfile(inputPath, '*.smrx')));
dataSmrxFile=fullfile(inputPath,smrxfilelist.name);
cedpath = 'F:\Sweyta\GRABS_Data\FinalPaperCodes\PreprocessingCode\Functions\spike2\CEDS64ML';

display(strcat('loading spike 2 file: ',dataSmrxFile));
if exist(fullfile(inputPath, 'smrx_signals.mat'),'file')
    load(fullfile(inputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx');
else
    [timing, channels_data] = process_spike2_two_cams_BS(cedpath, inputPath,dataSmrxFile, fsspike2,channels,minRunDuration,minSitDuration,ITITime,visStimDur,visStimITI,pupilSR);    
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
    save(fullfile(inputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx');
end

%% draw raw traces in visual parcell along with vis timestamps 
timeframes=1300:2000;
imaging_timestamps=timestamps.timaging(timeframes);
brainMapPlotTime=198.3; 
visual_times=timestamps.visStim(timestamps.visStim>imaging_timestamps(1) & timestamps.visStim<imaging_timestamps(end));
h1=figure; suptitle('TraceV1-GRABs');
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
linkaxes([ax1,ax2,ax3,ax4,ax5],'xy')
saveas(h1, fullfile(outputPath,'RawTracesDifferntRegressionMethods-GRABs')); 

h2=figure; suptitle('TraceV1-RCaMP');
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
linkaxes([ax1,ax2,ax3,ax4,ax5],'xy')

saveas(h2, fullfile(outputPath,'RawTracesDifferntRegressionMethods-RCaMP')); 

%% plot visual trial activity method 1
h1=figure; title('V1-GRABs-DFF');hold on; h2=figure; title('V1-RCaMP-DFF');hold on; V1Idx=2; 
h3=figure; title('V1-GRABs-ZScore');hold on; h4=figure; title('V1-RCaMP-Zscore');hold on;
V1data.blue=dFoF_parcells.blue(V1Idx,:);
V1data.green=dFoF_parcells.green(V1Idx,:);

%extract trials around visual stim 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%get difference df/f and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.blue-preEventMean).*100;%multiply by 100 to get percentages  
V1trial_aveDiffNorm.blue=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZNorm.blue=nanmean(ZNormTrial,2); 

preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:));
DiffNormtrial=(V1trial.green-preEventMean).*100; 
V1trial_aveDiffNorm.green=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZNorm.green=nanmean(ZNormTrial,2); 

set(0,'CurrentFigure',h1); ax1=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax2=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax11=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax12=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.green);box off; hold on; 


%% plot visual trial activity method 2
V1data.blue=dFoF_SVD_UV_parcells.blue(V1Idx,:);
V1data.green=dFoF_SVD_UV_parcells.green(V1Idx,:);

%extract trials around visual stim 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%get difference df/f and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.blue-preEventMean).*100;%multiply by 100 to get percentages  
V1trial_aveDiffNorm.blue=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZNorm.blue=nanmean(ZNormTrial,2); 

preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.green-preEventMean).*100; 
V1trial_aveDiffNorm.green=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZNorm.green=nanmean(ZNormTrial,2); 

set(0,'CurrentFigure',h1); ax3=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax4=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax13=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax14=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.green);box off; hold on; 

%% plot visual trial activity method 3
V1data.blue=dFoF_SVD_BS1_parcells.blue(V1Idx,:);
V1data.green=dFoF_SVD_BS1_parcells.green(V1Idx,:);

%extract trials around visual stim 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%get difference df/f and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.blue-preEventMean).*100;%multiply by 100 to get percentages  
V1trial_aveDiffNorm.blue=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZNorm.blue=nanmean(ZNormTrial,2); 

preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.green-preEventMean).*100; 
V1trial_aveDiffNorm.green=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZNorm.green=nanmean(ZNormTrial,2); 

set(0,'CurrentFigure',h1); ax5=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax6=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax15=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax16=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.green);box off; hold on; 

%% plot visual trial activity method 4
V1data.blue=dFoF_PixUV_parcells.blue(V1Idx,:);
V1data.green=dFoF_PixUV_parcells.green(V1Idx,:);

%extract trials around visual stim 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%get difference df/f and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.blue-preEventMean).*100;%multiply by 100 to get percentages  
V1trial_aveDiffNorm.blue=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZNorm.blue=nanmean(ZNormTrial,2); 

preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.green-preEventMean).*100; 
V1trial_aveDiffNorm.green=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZNorm.green=nanmean(ZNormTrial,2); 

set(0,'CurrentFigure',h1); ax7=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h2); ax8=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.green);box off; hold on; 
set(0,'CurrentFigure',h3); ax17=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.blue);box off; hold on; 
set(0,'CurrentFigure',h4); ax18=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.green);box off; hold on; 

%% plot visual trial activity method 5
V1data.blue=dFoF_PixBS1_parcells.blue(V1Idx,:);
V1data.green=dFoF_PixBS1_parcells.green(V1Idx,:);

%extract trials around visual stim 
eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
[V1trial.blue,timeStamps]= TrialTimeArrangeDff(V1data.blue,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
[V1trial.green,timeStamps]= TrialTimeArrangeDff(V1data.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);

%get difference df/f and average across trials 
preEventMean=nanmean(V1trial.blue(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.blue-preEventMean).*100;%multiply by 100 to get percentages  
V1trial_aveDiffNorm.blue=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.blue(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.blue-preEventMean)./preEventStd; 
V1trial_aveZNorm.blue=nanmean(ZNormTrial,2); 

preEventMean=nanmean(V1trial.green(1:preEventWin*fsimaging,:)); 
DiffNormtrial=(V1trial.green-preEventMean).*100; 
V1trial_aveDiffNorm.green=nanmean(DiffNormtrial,2); 
preEventStd=nanstd(V1trial.green(1:preEventWin*fsimaging,:)); 
ZNormTrial=(V1trial.green-preEventMean)./preEventStd; 
V1trial_aveZNorm.green=nanmean(ZNormTrial,2); 

set(0,'CurrentFigure',h1); ax9=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.blue);box off; hold on; 
legend([ax1,ax3,ax5,ax7,ax9],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1'}); legend boxoff; ylabel('Diff DFF'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h2); ax10=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveDiffNorm.green);box off; hold on; 
legend([ax2,ax4,ax6,ax8,ax10],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1'}); legend boxoff; ylabel('Diff DFF'); xlabel('Time (s)'); set(gca,'TickDir','out'); 

set(0,'CurrentFigure',h3); ax19=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.blue);box off; hold on; 
legend([ax11,ax13,ax15,ax17,ax19],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1'}); legend boxoff; ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 
set(0,'CurrentFigure',h4); ax20=plot((-preEventWin:1/fsimaging:(postEventWin-1/fsimaging))',V1trial_aveZNorm.green);box off; hold on; 
legend([ax12,ax14,ax16,ax18,ax20],{'Uncorrected','SpatialSVD-UV','SpatialSVD-BS1','PixelwiseRegress-UV','PixelwiseRegress-BS1'}); legend boxoff; ylabel('ZScore'); xlabel('Time (s)'); set(gca,'TickDir','out'); 

saveas(h1,fullfile(outputPath,'VisualStimAverage-V1-GRABs-DFF')); saveas(h1,fullfile(outputPath,'VisualStimAverage-V1-GRABs-DFF.pdf')); 
saveas(h2,fullfile(outputPath,'VisualStimAverage-V1-RCaMP-DFF')); saveas(h2,fullfile(outputPath,'VisualStimAverage-V1-RCaMP-DFF.pdf')); 
saveas(h3,fullfile(outputPath,'VisualStimAverage-V1-GRABs-ZScore')); saveas(h1,fullfile(outputPath,'VisualStimAverage-V1-GRABs-ZScore.pdf')); 
saveas(h4,fullfile(outputPath,'VisualStimAverage-V1-RCaMP-ZScore')); saveas(h2,fullfile(outputPath,'VisualStimAverage-V1-RCaMP-ZScore.pdf')); 
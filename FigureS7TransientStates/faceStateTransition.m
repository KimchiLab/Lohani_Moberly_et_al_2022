%% Figure S7 This script calculates event triggered averages at face state transitions (low face to high face during sustained quiescence)
%Sweyta Lohani 2020 
%% inputs and outputs 
inputFolder='W:\GRABS_Data\Analyzed_SVDMethodPatch14\';% where input data are located 
outputFolder='F:\Figures\TransitionStates\Face';
MiceAnalyze=[{'DualMice\grabAM10\imaging with 575 excitation'},{'DualMice\grabAM09\imaging with 575 excitation'},{'DualMice\grabAM08\imaging with 575 excitation'},{'DualMice\grabAM07\imaging with 575 excitation'},{'DualMice\grabAM06\imaging with 575 excitation'},{'DualMice\grabAM05\imaging with 575 excitation'}];
Condition='NoDrug'; %'NoDrug','PreDrug','PostDrug'
%% user-selected input parameters
params.fsspike2=5000;%spike2 sampling frequency
params.fsimaging=10;%imaging sampling frequency
params.fspupilcam=10;%pupil/face cam sampling frequency
params.signalsExtraction.sigs = 'RCaMP_AC'; % 'blueuv' or 'RCaMP_AC'

%analysis window for events (face and pupil) during quiescence
params.TimeSinceSitOn=10;%for sustained quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=10;%for sustained quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=10;%for any state, minimum time since any event onset/offset such as airpuff
params.minSitDuration=5;%minimum sit duration during sustained quiescnece state
params.minArousalDuration=0; %minimum face/pupil arousal state (high or low arousal)
params.TimeSinceArousalOn=0;%for face/pupil arousal state, minimum time since  onset
params.TimeBeforeArousalOff=0;%for face/pupil arousal state, minimum time before  offset
params.TimeBeforeArousalOn=4;%for face/pupil arousal state, minimum time of low arousal before high arousal state

params.preEventWin=4;%pre event window for face and pupil trials in seconds
params.postEventWin=4;%post-event window for face and pupil trials in seconds
params.baselineWin=[-4 -2]; %use this period to calculate baseline during within trial normalization for face and pupil trials
params.statWin=1; %calculate summary stats such as mean or max in this moving window in seconds
params.statTime=0; %calculate statistical differences and significance at this time window post event start
params.wheelDownsSampleFactor=100;%down sample wheel by this rate

%EEG parameters
params.win=1; params.step=0.05;%window and step size
params.Fs=5000; % sampling frequency
params.fpass=[1 100]; % frequencies of interest
params.tapers=[2 3]; % tapers
params.trialave=0; % average over trials
%BrainRegions
params.Regions=[{'LeftV1'}; {'LeftS1b'};{'LeftM2'}];
params.RegionColor=[{[0.9 0.19 0.2]},{[1 0.55 0.1]},{[0.82 0 1]}];

%example mouse and session
params.Mouse=1; params.Session=3; %example mouse and session for representative behavior figures 
%% add path for functions
addpath(genpath('./Functions'));
if ~exist(outputFolder,'dir'),mkdir(outputFolder); end
%% get parcell idx for both hemishpers
load('parcells_updated121519.mat'); parcells=parcells_new;
Idx.visual=1:1:16; %visual parcells
Idx.PosteriorParietal=[31,32];% PPC
Idx.RSL_AgL=[17,18,19,20];%retrosplenial cortex
Idx.Somat=33:1:48;%somatosensory areas
Idx.FrontalMotor=49:1:52;%frontal-moto area
Idx.auditory=[27,28,29,30];%auditory areas
CombinedParcellIdx=[Idx.visual,Idx.RSL_AgL,Idx.auditory,Idx.PosteriorParietal,Idx.Somat,Idx.FrontalMotor];
parcellnames=parcells.names(CombinedParcellIdx);
V1Idx=2; S1bIdx=28; M2Idx=46; %updated parcell idx for left V1, S1, and M2

%% combine data across mice
if ~exist(fullfile(outputFolder,'SpontCombinedVars.mat'),'file')
    [Comb,DataS]= combineData_spontFacePupil(inputFolder,outputFolder,MiceAnalyze,params,Condition);
    save(fullfile(outputFolder,'SpontCombinedVars'),'Comb','-v7.3');
    save(fullfile(outputFolder,'SpontSessionVars'),'DataS','-v7.3');
else
    load(fullfile(outputFolder,'SpontCombinedVars'),'Comb');
    load(fullfile(outputFolder,'SpontSessionVars'),'DataS');
end

%% make population averages for face state transition 
Blue=[]; Green=[];
%blue imaging face 
h1=figure; h2=figure;color='Blue';  eventName='face';NormType='DiffNorm'; ylab='DFF';ylimit=[-0.5 1];
for i=1:length(params.Regions)
    currRegion=params.Regions{i};
    [Blue]=makePopAverages(Blue,Comb.blue, eventName,NormType,currRegion,params);
    [h1,h2,handle(i),handle1(i),Blue]=plotPopAverage(h1,h2,Blue,eventName,NormType,currRegion,params.fsimaging,params.preEventWin, params.postEventWin,ylab,ylimit,color,params.RegionColor{i});
end
set(0, 'CurrentFigure', h1);line([0,0],ylim,'linestyle','--','Color','r');legend([handle(1:end)],params.Regions); legend boxoff; hold off
saveas (h1,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-Animals')));saveas (h1,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-Animals.pdf')));
set(0, 'CurrentFigure', h2);line([0,0],ylim,'linestyle','--','Color','r');legend([handle1(1:end)],params.Regions); legend boxoff; hold off
saveas (h2,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-ExampleMouse-Trials')));saveas (h2,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-ExampleMouse-Trials.pdf')));
%do stats on the max response at specified window from face onset
frame=((params.preEventWin+params.statTime)*params.statWin+1);
peakStats=[];
for i=1:length(params.Regions)
    currRegion=params.Regions{i};
    peakStats=cat(2,peakStats,Blue.face.DiffNorm.(strcat(currRegion,'_MaxWin'))(:,frame));
end
Blue.face.PeakStats=peakStats;
h3=figure;plot(Blue.face.PeakStats','-ko','MarkerFaceColor','k');ylim([0 3]);ylabel(ylab); xticks([1  2 3]); xticklabels(params.Regions); title(strcat('PairedPlot',eventName,color)); box off;hold on;
for i=1:length(params.Regions)
    currMedian=median(Blue.face.PeakStats(:,i));
    hLine=plot(i,currMedian,'+'); hLine.MarkerEdgeColor = params.RegionColor{i}; %draw horizontal bars at median
end
hold off;
saveas (h3,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-PairedDataAnimals')));saveas (h3,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-PairedDataAnimals.pdf')));

%get data for all parcells
ylimBrainMap=[-1.5 1.5];
[figure5,figure6, Blue]=makePopAverages_allparcells(Blue,Comb.blue, eventName,NormType,params,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,parcellnames,ylimBrainMap);
set(0, 'CurrentFigure', figure5);title(strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF'));
saveas (figure5,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF-BrainMap'))); saveas (figure5,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF-BrainMap.pdf')));
set(0, 'CurrentFigure', figure6);title(strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF'));
saveas (figure6,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-Animal-PeakDFF-BarGraph'))); saveas (figure6,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-Animal-PeakDFF-BarGraph.pdf')));
%% green imaging face
if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
    h1=figure; h2=figure; color='Green';  eventName='face';NormType='DiffNorm'; ylab='DFF';ylimit=[-0.5 2];
    for i=1:length(params.Regions)
        currRegion=params.Regions{i};
        [Green]=makePopAverages(Green,Comb.green, eventName,NormType,currRegion,params);
        [h1,h2,handle(i),handle1(i),Green]=plotPopAverage(h1,h2,Green,eventName,NormType,currRegion,params.fsimaging,params.preEventWin, params.postEventWin,ylab,ylimit,color,params.RegionColor{i});
    end
    set(0, 'CurrentFigure', h1);line([0,0],ylim,'linestyle','--','Color','r');legend([handle(1:end)],params.Regions); legend boxoff; hold off
    saveas (h1,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-Animals')));saveas (h1,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-Animals.pdf')));
    set(0, 'CurrentFigure', h2);line([0,0],ylim,'linestyle','--','Color','r');legend([handle1(1:end)],params.Regions); legend boxoff; hold off
    saveas (h2,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-ExampleMouse-Trials')));saveas (h2,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-ExampleMouse-Trials.pdf')));
    
    %do stats on the max response at specified window from face onset
    frame=((params.preEventWin+params.statTime)*params.statWin+1);
    peakStats=[];
    for i=1:length(params.Regions)
        currRegion=params.Regions{i};
        peakStats=cat(2,peakStats,Green.face.DiffNorm.(strcat(currRegion,'_MaxWin'))(:,frame));
    end
    Green.face.PeakStats=peakStats;
    h3=figure;plot(Green.face.PeakStats','-ko','MarkerFaceColor','k');ylim([0 3]); ylabel(ylab); xticks([1  2 3]); xticklabels(params.Regions); title(strcat('PairedPlot',eventName,color));box off;hold on;
    for i=1:length(params.Regions)
        currMedian=median(Green.face.PeakStats(:,i));
        hLine=plot(i,currMedian,'+'); hLine.MarkerEdgeColor = params.RegionColor{i};%draw horizontal bars at median
    end
    hold off;
    saveas (h3,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-PairedDataAnimals')));saveas (h3,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAverage-PairedDataAnimals.pdf')));
    %get data for all parcells
    ylimBrainMap=[-2.5 2.5];
    [figure7,figure8, Green]=makePopAverages_allparcells(Green,Comb.green, eventName,NormType,params,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,parcellnames,ylimBrainMap);
    set(0, 'CurrentFigure', figure7);title(strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF'));
    saveas (figure7,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF-BrainMap'))); saveas (figure7,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF-BrainMap.pdf')));
    set(0, 'CurrentFigure', figure8);title(strcat(color,'-',eventName,'-','DiffNorm','-PopAve-Animal-PeakDFF'));
    saveas (figure8,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-Animal-PeakDFF-BarGraph'))); saveas (figure8,fullfile(outputFolder,strcat(color,'-',eventName,'-','DiffNorm','-Animal-PeakDFF-BarGraph.pdf')));
end
save(fullfile(outputFolder,'summaryData.mat'),'Blue','Green'); %save all the summary data including stats

%% plot behavioral results
% plot face run speed, pupil, face and EEG
eventName='face'; preEventWin=params.preEventWin;
[h5,h6]=plotBehaviorAverage(Comb,eventName,preEventWin,params);
saveas (h5,fullfile(outputFolder,strcat('Behavior-',eventName,'-AnimalAve')));saveas (h5,fullfile(outputFolder,strcat('Behavior-',eventName,'-AnimalAve.pdf')));
saveas (h6,fullfile(outputFolder,strcat('Behavior-',eventName,'-ExampleMouseTrialAve')));saveas (h6,fullfile(outputFolder,strcat('Behavior-',eventName,'-ExampleMouseTrialAve.pdf')));

%% This script makes example state figures in Figure 3a
%written by Sweyta Lohani 2020
close all; clear all
%% user defined folder inputs
figuresFolder='F:\Figures';%Where figures will be saved 
inputFolder='W:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice';%input data moler
MiceAnalyze=[{'grabAM05\imaging with 575 excitation\'}];%mouse to analyze 
Condition='NoDrug'; 
%% user-selected input parameters
params.signalsExtraction= 'RCaMP_AC'; % 'blueuv' or 'RCaMP_AC'
params.fsimaging=10;%imaging sampling rate
params.fspupilcam=10; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=10;%for any state, minimum time since any event onset/offset
params.minRunDuration=5;% minimum run duration during locomotion state
params.minArousalDuration=5; %minimum face/pupil arousal state (high or low arousal)
params.TimeBeforeArousalOn=4;%for face/pupil arousal state, minimum time of no/low arousal before high arousal state
params.minSitDuration=5;%minimum sit duration during quiescnece state
params.win=3; params.step=0.5; %window and steps for EEG 

%% add functions
addpath(genpath('./Functions'));
if ~exist(figuresFolder),mkdir(figuresFolder); end

%% get parcell indices for parcells we care about in the left and right hemisphere
load('parcells_updated121519.mat'); parcells=parcells_new;
Idx.visual=2:2:16; %visual parcells
Idx.PosteriorParietal=32;% PPC
Idx.RSL_AgL=[18,20];%retrosplenial cortex
Idx.Somat=34:2:48;%somatosensory areas
Idx.FrontalMotor=50:2:52;%frontal-moto area
Idx.auditory=[28,30];%auditory areas
CombinedParcellIdx=[Idx.visual,Idx.RSL_AgL,Idx.auditory,Idx.PosteriorParietal,Idx.Somat,Idx.FrontalMotor];

%get the size of each parcell so we can scale each cell in the color map accrdingly
row=size(parcells.indicators,1); col=size(parcells.indicators,2);
pixelSizeParcell=sum(reshape(parcells.indicators(:,:,CombinedParcellIdx),row*col,length(CombinedParcellIdx)),1); %get the size of each parcell
totalBrainPixels=sum(pixelSizeParcell);
pixelpropSizeParcell=pixelSizeParcell/totalBrainPixels;
parcellnames=parcells.names(CombinedParcellIdx);
V1Idx=1; S1bIdx=14; M2Idx=23; %updated parcell idx for left V1, S1, and M2

%% get states for one example session and folder
animal=1; folder =1;
mainDir1 =fullfile(inputFolder,MiceAnalyze{animal});
folders=dir(mainDir1);
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
DirFolders= folders(dirFlags);
noDrugFlag=~contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
preDrugFlag=contains({DirFolders.name},'preDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
postDrugFlag=contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
switch Condition
    case 'NoDrug'
        DirFolders=DirFolders(noDrugFlag);
    case 'PreDrug'
        DirFolders=DirFolders(preDrugFlag);
    case 'PostDrug'
        DirFolders=DirFolders(postDrugFlag);
end

currfolder=fullfile(mainDir1, DirFolders(folder).name);

% load spike2 and pupil/face data
load(fullfile(currfolder,'smrx_signals.mat'),'channels_data','timestamps','timing')
files=dir(currfolder);
for t=1:length(files)
    if contains(files(t).name,'proc','IgnoreCase',true)
        files=files(t);break;
    end
end
load(fullfile(currfolder,files.name),'proc')
pupil_Norm=proc.output.pupilNorm;
face_Norm=proc.output.facePC1CorrNorm;
wheel_speed = channels_data.wheelspeed(1:end-1);%last element is NaN
EEG=channels_data.EEG(1:end-1);%last element is NaN

% get pupil wheel and imaging times
wheel_time = (1:length(wheel_speed))/params.fsspike2;
EEG_time=(1:length(EEG))/params.fsspike2;
imaging_time = timestamps.timaging;
pupil_time = timing.pupilcamstart(1:length(pupil_Norm));

%% make output figures folder
indivFigureFolder=fullfile(figuresFolder,MiceAnalyze{animal},DirFolders(folder).name);
if ~exist(indivFigureFolder,'dir'),mkdir(indivFigureFolder), end
%% calculate EEG bandpower
movingwin=[params.win params.step]; % set the moving window params
pms.Fs=5000; % sampling frequency
pms.fpass=[1 100]; % frequencies of interest
pms.tapers=[5 9]; % tapers
[S1,t,f]=mtspecgramc(EEG,movingwin,pms);
baseS1=mean(S1,1); stdS1=std(S1,0,1);  %normalize to mean activity
S1norm=(S1-baseS1)./stdS1;
lowPower=mean(S1norm(:,f<10),2);
highPower=mean(S1norm(:,f>30),2);
EEG_SpecTime=t; 
freq=f; 

%get amplitude envelope of EEG in 1-10 Hz 
[filteredEEG] = eegfilt(EEG',params.fsspike2,1,10); 
EEG_env=abs(hilbert(filteredEEG));

%% get locomotion on/off and quiescence on off times
%locomotion periods should be at least some criterion s long with some criterion s since locomotion
%onset, some criterion s before locomotion offset, excluding any events (airpuff/stim)
wheelOn=timing.allwheelon;
wheelOff=timing.allwheeloff;

%find wheel on and wheel off times during imaging period only
minRunDur=params.minRunDuration+params.TimeSinceLocOn+params.TimeBeforeLocOff; %minimum actual locomotion duration including time since locomotion onset, time before locomotion offset and the minimum time period for data analysis
idx=wheelOn<(imaging_time(end)) & wheelOff>(imaging_time(1));
wheelOn_t1=wheelOn(idx);
wheelOff_t1=wheelOff(idx);

for whe=1:length(wheelOn_t1)
    if wheelOn_t1(whe)<imaging_time(end) && wheelOff_t1(whe)>imaging_time(end)
        wheelOff_t1(whe)=imaging_time(end)-1;%if locomotion starts before end of imaging but continues after, only extract state until imaging time end minus a second
    end
    
    if wheelOn_t1(whe)<imaging_time(1) && wheelOff_t1(whe)>imaging_time(1)
        wheelOn_t1(whe)=imaging_time(1)+1;%if locomotion starts before start of imaging but continues after, only extract state from imaging time start plus a second
    end
end

wheelOn_t1=(wheelOn_t1(:))'; wheelOff_t1=(wheelOff_t1(:))';

%find wheel on off  times when airpuffs are not given
if ~isempty(timing.airpuffstart)&& ~isempty(wheelOn_t1)
    allEvts=sort(timing.airpuffstart,'ascend');
    allEvts=allEvts(:)';
    allEvtsPre=allEvts-params.TimeSinceEvent;
    allEvtsPost=allEvts+params.TimeSinceEvent;
    
    newwheelOn=sort([wheelOn_t1,allEvtsPre,allEvtsPost],'ascend');
    newwheelOff=sort([wheelOff_t1,allEvtsPre,allEvtsPost],'ascend');
    
    Index=nan(1,length(newwheelOn));
    for r=1:length(newwheelOn)
        tmp1=newwheelOn(r);
        tmp2=newwheelOff(r);
        
        tmp3= sum(tmp1>=allEvtsPre & tmp2<=allEvtsPost);
        if sum(tmp3)==0
            Index(r)=r;
        end
    end
    wheelOn_int=newwheelOn(Index(~isnan(Index)));
    wheelOff_int=newwheelOff(Index(~isnan(Index)));
else
    wheelOn_int=wheelOn_t1;
    wheelOff_int=wheelOff_t1;
end

%makes sure the state is at least as long as the minimum run duration
idx1=find((wheelOff_int-wheelOn_int)>=(minRunDur));
wheelOn_int1=wheelOn_int(idx1);
wheelOff_int1=wheelOff_int(idx1);

%finalize the times to get sustained state only
wheelOn_final=wheelOn_int1+params.TimeSinceLocOn;
wheelOff_final=wheelOff_int1-params.TimeBeforeLocOff;

%% get wheel on/off times for state transition points only (timestamps used for transition average figures)
wheelOn_tran=timing.wheelOn; 
wheelOff_tran=timing.wheelOff; 

%% queiscence should be at least some criterion s long with some criterion s since locomotion offset
%and some criterion s before subsequent locomotion onset,excluding any events (airpuff/stim)
sitOn=[0;wheelOff(1:end-1)]; %use 0 as the first sit on time;
sitOff=wheelOn;%use wheelOn times as sit off times;

%find sit on and sit off times during imaging period only
minSitDur=params.minSitDuration+params.TimeSinceSitOn+params.TimeBeforeSitOff; %actual minimum sit duration accouting for the onset time, offset time and minimum duration of the sustained quiescence epoch  used for analysis

idx=sitOn<(imaging_time(end)) & sitOff>(imaging_time(1));
sitOn_t1=sitOn(idx);
sitOff_t1=sitOff(idx);

for whe=1:length(sitOn_t1)
    if sitOn_t1(whe)<imaging_time(end) && sitOff_t1(whe)>imaging_time(end)
        sitOff_t1(whe)=imaging_time(end)-1;%if queiscence starts before end of imaging but continues after, only extract state until imaging time end minus a second
    end
    
    if sitOn_t1(whe)<imaging_time(1) && sitOff_t1(whe)>imaging_time(1)
        sitOn_t1(whe)=imaging_time(1)+1;%if queiscence starts before start of imaging but continues after, only extract state from imaging time start plus a second
    end
end

%remove any quiescence period where mouse's speed is above 0.03 m/s
% because sometimes these epcohs,especially short ones, get missed by locomotion changepoint algorithm
speedThres=0.03;
tmpOn=sitOn_t1; tmpOff=sitOff_t1;
for rr=1:length(sitOn_t1)
    highSpeedIdx=find(abs(wheel_speed)>speedThres & wheel_time>sitOn_t1(rr) & wheel_time<sitOff_t1(rr));
    if ~isempty(highSpeedIdx)
        firstIdx=wheel_time(highSpeedIdx(1))-0.1; lastIdx=wheel_time(highSpeedIdx(end))+0.1;
        tmpOff(rr)=firstIdx; tmpOff(end+1)=sitOff_t1(rr); tmpOn(end+1)=lastIdx;
    end
end
tmpOn=sort(tmpOn,'ascend'); tmpOff=sort(tmpOff,'ascend');

sitOn_t1=(tmpOn(:))'; sitOff_t1=(tmpOff(:))';

%find sit on off  times when airpuffs are not given
allEvts=sort(timing.airpuffstart,'ascend');
allEvts=allEvts(:)';
if ~isempty(allEvts)&& ~isempty(sitOn_t1)
    allEvtsPre=allEvts-params.TimeSinceEvent;
    allEvtsPost=allEvts+params.TimeSinceEvent;
    
    newSitOn=sort([sitOn_t1,allEvtsPre,allEvtsPost],'ascend');
    newSitOff=sort([sitOff_t1,allEvtsPre,allEvtsPost],'ascend');
    
    Index=nan(1,length(newSitOn));
    for r=1:length(newSitOn)
        tmp1=newSitOn(r);
        tmp2=newSitOff(r);
        
        tmp3= sum(tmp1>=allEvtsPre & tmp2<=allEvtsPost);
        if sum(tmp3)==0
            Index(r)=r;
        end
    end
    sitOn_int=newSitOn(Index(~isnan(Index)));
    sitOff_int=newSitOff(Index(~isnan(Index)));
else
    sitOn_int=sitOn_t1;
    sitOff_int=sitOff_t1;
end

%makes sure state is at least as long as the minimum sit duration
idx1=find((sitOff_int-sitOn_int)>=(minSitDur));
sitOn_int1=sitOn_int(idx1);
sitOff_int1=sitOff_int(idx1);

%finalize the times to get sustained state only
sitOn_final=sitOn_int1+params.TimeSinceSitOn;
sitOff_final=sitOff_int1-params.TimeBeforeSitOff;

%% do change point detection on pupil and face to get pupil high/low arousal or face high/low movement times during sustained quiescence state
%get Z-thresholds based on pupil data during quiescence, when mouse
%isn't moving and when aripuffs are not given
b1DoPlot=0; blDoPlotDuration=1:4000; smoothWin=1;
pupilTime_Idx=cell(1,length(sitOn_int));
for st=1:length(sitOn_int)
    pupilTime_Idx{st}=find(pupil_time>sitOn_int(st) & pupil_time <sitOff_int(st));
end
pupilTime_quiescence=cell2mat(pupilTime_Idx');
pupil_qui_times=pupil_time(pupilTime_quiescence);
pupil_quiescence=pupil_Norm(pupilTime_quiescence);
zthres_High=quantile(pupil_quiescence,0.60);
zthres_Low=quantile(pupil_quiescence,0.40);

%get on and off timestamps for high and low arousal based on pupil
[~,Pupil_HighArousal_OnTStamp,Pupil_HighArousal_OffTStamp ] =changepoints(pupil_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
[~,Pupil_LowArousal_OnTStamp,Pupil_LowArousal_OffTStamp ] =changepoints(-pupil_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
%get z thresholds based on face data during quiescence only
b1DoPlot=0; blDoPlotDuration=1:4000; smoothWin=1;
face_quiescence=face_Norm(pupilTime_quiescence);
zthres_High=quantile(face_quiescence,0.60);
zthres_Low=quantile(face_quiescence,0.40);

%get on and off timestamps for high and low face movment
[~,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
[~,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);

% determine that pupil and face high/low arousal/movement time are at least minimum criterion seconds long
idx1=find((Pupil_HighArousal_OffTStamp-Pupil_HighArousal_OnTStamp)>=params.minArousalDuration);
Pupil_HighArousal_OnT_int=Pupil_HighArousal_OnTStamp(idx1); Pupil_HighArousal_OffT_int=Pupil_HighArousal_OffTStamp(idx1);
idx2=find((Pupil_LowArousal_OffTStamp-Pupil_LowArousal_OnTStamp)>=params.minArousalDuration);
Pupil_LowArousal_OnT_int=Pupil_LowArousal_OnTStamp(idx2); Pupil_LowArousal_OffT_int=Pupil_LowArousal_OffTStamp(idx2);
idx3=find((Face_HighArousal_OffTStamp-Face_HighArousal_OnTStamp)>=params.minArousalDuration);
Face_HighArousal_OnT_int=Face_HighArousal_OnTStamp(idx3); Face_HighArousal_OffT_int=Face_HighArousal_OffTStamp(idx3);
idx4=find((Face_LowArousal_OffTStamp-Face_LowArousal_OnTStamp)>=params.minArousalDuration);
Face_LowArousal_OnT_int=Face_LowArousal_OnTStamp(idx4); Face_LowArousal_OffT_int=Face_LowArousal_OffTStamp(idx4);

% get pupil on and face on times if both on and off times occur entirely during sustained queiscence states identified in the previous step
toDelete=ones(1,length(Pupil_HighArousal_OnT_int));
for rj=1:length(Pupil_HighArousal_OnT_int)
    tmp = find (Pupil_HighArousal_OnT_int(rj)>=sitOn_final & Pupil_HighArousal_OffT_int(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Pupil_HighArousal_On_final=Pupil_HighArousal_OnT_int(~toDelete);
% Pupil_HighArousal_Off_final=Pupil_HighArousal_OffT_int(~toDelete);
Pupil_HighArousal_On_final=Pupil_HighArousal_OnT_int;
Pupil_HighArousal_Off_final=Pupil_HighArousal_OffT_int;

toDelete=ones(1,length(Pupil_LowArousal_OnT_int));
for rj=1:length(Pupil_LowArousal_OnT_int)
    tmp = find (Pupil_LowArousal_OnT_int(rj)>=sitOn_final & Pupil_LowArousal_OffT_int(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Pupil_LowArousal_On_final=Pupil_LowArousal_OnT_int(~toDelete);
% Pupil_LowArousal_Off_final=Pupil_LowArousal_OffT_int(~toDelete);
Pupil_LowArousal_On_final=Pupil_LowArousal_OnT_int;
Pupil_LowArousal_Off_final=Pupil_LowArousal_OffT_int;

toDelete=ones(1,length(Face_HighArousal_OnT_int));
for rj=1:length(Face_HighArousal_OnT_int)
    tmp = find (Face_HighArousal_OnT_int(rj)>=sitOn_final & Face_HighArousal_OffT_int(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Face_HighArousal_On_final=Face_HighArousal_OnT_int(~toDelete);
% Face_HighArousal_Off_final=Face_HighArousal_OffT_int(~toDelete);
Face_HighArousal_On_final=Face_HighArousal_OnT_int;
Face_HighArousal_Off_final=Face_HighArousal_OffT_int;

toDelete=ones(1,length(Face_LowArousal_OnT_int));
for rj=1:length(Face_LowArousal_OnT_int)
    tmp = find (Face_LowArousal_OnT_int(rj)>=sitOn_final & Face_LowArousal_OffT_int(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Face_LowArousal_On_final=Face_LowArousal_OnT_int(~toDelete);
% Face_LowArousal_Off_final=Face_LowArousal_OffT_int(~toDelete);
Face_LowArousal_On_final=Face_LowArousal_OnT_int;
Face_LowArousal_Off_final=Face_LowArousal_OffT_int;


%% get pupil/face on/off times for state transition points only (timestamps used for transition average figures)

% get pupil on and face on times if both on and off times occur entirely during sustained queiscence states identified in the previous step
toDelete=ones(1,length(Pupil_HighArousal_OnTStamp));
for rj=1:length(Pupil_HighArousal_OnTStamp)
    tmp = find (Pupil_HighArousal_OnTStamp(rj)>=sitOn_final & Pupil_HighArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Pupil_HighArousal_On_int=Pupil_HighArousal_OnTStamp(~toDelete);
% Pupil_HighArousal_Off_int=Pupil_HighArousal_OffTStamp(~toDelete);
Pupil_HighArousal_On_int=Pupil_HighArousal_OnTStamp;
Pupil_HighArousal_Off_int=Pupil_HighArousal_OffTStamp;

toDelete=ones(1,length(Pupil_LowArousal_OnTStamp));
for rj=1:length(Pupil_LowArousal_OnTStamp)
    tmp = find (Pupil_LowArousal_OnTStamp(rj)>=sitOn_final & Pupil_LowArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Pupil_LowArousal_On_int=Pupil_LowArousal_OnTStamp(~toDelete);
% Pupil_LowArousal_Off_int=Pupil_LowArousal_OffTStamp(~toDelete);
Pupil_LowArousal_On_int=Pupil_LowArousal_OnTStamp;
Pupil_LowArousal_Off_int=Pupil_LowArousal_OffTStamp;

toDelete=ones(1,length(Face_HighArousal_OnTStamp));
for rj=1:length(Face_HighArousal_OnTStamp)
    tmp = find (Face_HighArousal_OnTStamp(rj)>=sitOn_final & Face_HighArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Face_HighArousal_On_int=Face_HighArousal_OnTStamp(~toDelete);
% Face_HighArousal_Off_int=Face_HighArousal_OffTStamp(~toDelete);
Face_HighArousal_On_int=Face_HighArousal_OnTStamp;
Face_HighArousal_Off_int=Face_HighArousal_OffTStamp;

toDelete=ones(1,length(Face_LowArousal_OnTStamp));
for rj=1:length(Face_LowArousal_OnTStamp)
    tmp = find (Face_LowArousal_OnTStamp(rj)>=sitOn_final & Face_LowArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
% Face_LowArousal_On_int=Face_LowArousal_OnTStamp(~toDelete);
% Face_LowArousal_Off_int=Face_LowArousal_OffTStamp(~toDelete);
Face_LowArousal_On_int=Face_LowArousal_OnTStamp;
Face_LowArousal_Off_int=Face_LowArousal_OffTStamp;

%get high face/pupil states that are preceded by low arousal state
finalIndex=zeros(1,length(Face_HighArousal_On_int));
for ll=1:length(Face_HighArousal_On_int)
    currOnPoint=Face_HighArousal_On_int(ll);
    tmp=currOnPoint-Face_LowArousal_Off_int;
    idxtmp=find(tmp>=0 &tmp<=0.5); %find states where low arousal ends at less than 0.5s before high arousal onset;
    lowarousallength=Face_LowArousal_Off_int(idxtmp)-Face_LowArousal_On_int(idxtmp);
    if length(lowarousallength)>1 || isempty(lowarousallength), continue, end
    if lowarousallength>=params.TimeBeforeArousalOn
        finalIndex(ll)=1;
        %Face_HighArousal_On_int(ll)=Face_LowArousal_Off_int(idxtmp); %replace aorusal on time with low arousal off time plus two samples
    end
end
Face_HighArousal_On_tran=Face_HighArousal_On_int(logical(finalIndex));
Face_HighArousal_Off_tran=Face_HighArousal_Off_int(logical(finalIndex));

finalIndex=zeros(1,length(Pupil_HighArousal_On_int));
for ll=1:length(Pupil_HighArousal_On_int)
    currOnPoint=Pupil_HighArousal_On_int(ll);
    tmp=currOnPoint-Pupil_LowArousal_Off_int;
    idxtmp=find(tmp>=0 &tmp<=0.5); %find states where low arousal ends at less than 0.5s before high arousal onset;
    lowarousallength=Pupil_LowArousal_Off_int(idxtmp)-Pupil_LowArousal_On_int(idxtmp);
    if length(lowarousallength)>1 || isempty(lowarousallength), continue, end
    if lowarousallength>=params.TimeBeforeArousalOn
        finalIndex(ll)=1;
        Pupil_HighArousal_On_int(ll)=Pupil_LowArousal_Off_int(idxtmp); %replace aorusal on time with low arousal off time
    end
end
Pupil_HighArousal_On_tran=Pupil_HighArousal_On_int(logical(finalIndex));
Pupil_HighArousal_Off_tran=Pupil_HighArousal_Off_int(logical(finalIndex));


%% plot on and off times for transient and sustained states as well as EEG 
%plot locomotion state on off times
h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified states
ax1=subplot(6,1,1);plot(wheel_time,wheel_speed,'k');ylimits=ylim; hold on; title('LocomotionStates');
for tt=1:length(wheelOn_final)    
    rectangle('Position',[wheelOn_final(tt) ylimits(1) wheelOff_final(tt)-wheelOn_final(tt) ylimits(2)-ylimits(1)],'FaceColor',[0 0 1 0.2],'EdgeColor','none')%locomotion sustained state 
end
for tt=1:length(wheelOn_tran) %lcoomotion transition points 
plot([wheelOn_tran(tt),wheelOn_tran(tt)], ylimits, 'g');
plot([wheelOff_tran(tt),wheelOff_tran(tt)], ylimits, 'r');
end 
box off; set(gca,'TickDir','out');

for tt=1:length(sitOn_final)
    rectangle('Position',[sitOn_final(tt) ylimits(1) sitOff_final(tt)-sitOn_final(tt) ylimits(2)-ylimits(1)],'FaceColor',[1 0 0 0.2],'EdgeColor','none'); %queiesnce sustained state
end

%plot pupil on off times
ax2=subplot(6,1,2);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilStates');
for tt=1:length(Pupil_HighArousal_On_final)
    rectangle('Position',[Pupil_HighArousal_On_final(tt) ylimits(1) Pupil_HighArousal_Off_final(tt)-Pupil_HighArousal_On_final(tt) ylimits(2)-ylimits(1)],'FaceColor',[0 0 1 0.2],'EdgeColor','none')%pupil high sustained state 
end

for tt=1:length(Pupil_HighArousal_On_tran)%pupil transition points 
    plot([Pupil_HighArousal_On_tran(tt),Pupil_HighArousal_On_tran(tt)], ylimits, 'g');
    plot([Pupil_HighArousal_Off_tran(tt),Pupil_HighArousal_Off_tran(tt)], ylimits, 'r');
end

for tt=1:length(Pupil_LowArousal_On_final)
     rectangle('Position',[Pupil_LowArousal_On_final(tt) ylimits(1) Pupil_LowArousal_Off_final(tt)-Pupil_LowArousal_On_final(tt) ylimits(2)-ylimits(1)],'FaceColor',[1 0 0 0.2],'EdgeColor','none')%pupil low sustained state 
end
box off; set(gca,'TickDir','out');

%plot face on off times 
ax3=subplot(6,1,3);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceStates');
for tt=1:length(Face_HighArousal_On_final)
    rectangle('Position',[Face_HighArousal_On_final(tt) ylimits(1) Face_HighArousal_Off_final(tt)-Face_HighArousal_On_final(tt) ylimits(2)-ylimits(1)],'FaceColor',[0 0 1 0.2],'EdgeColor','none')%pupil high sustained state 
end

for tt=1:length(Face_HighArousal_On_tran)%pupil transition points 
    plot([Face_HighArousal_On_tran(tt),Face_HighArousal_On_tran(tt)], ylimits, 'g');
    plot([Face_HighArousal_Off_tran(tt),Face_HighArousal_Off_tran(tt)], ylimits, 'r');
end

for tt=1:length(Face_LowArousal_On_final)
     rectangle('Position',[Face_LowArousal_On_final(tt) ylimits(1) Face_LowArousal_Off_final(tt)-Face_LowArousal_On_final(tt) ylimits(2)-ylimits(1)],'FaceColor',[1 0 0 0.2],'EdgeColor','none')%pupil low sustained state 
end
box off; set(gca,'TickDir','out');

% EEG spectrogram
ax4=subplot(6,1,4);imagesc(EEG_SpecTime,freq,S1norm'); axis xy; caxis([-1 1]); title('EEG'); 
box off; set(gca,'TickDir','out');

linkaxes([ax1, ax2, ax3,ax4],'x');

saveas(h,fullfile(indivFigureFolder,'BehavioralStatesExample.fig')); 

xlim([1250 1675])
%% This script calculate correlations (cross-correlation with specified time lag) between parcells for different color combinations (blue-blue, green-green)
%or within parcell for two colors (Blue and green). All time are in seconds
%written by Sweyta Lohani 2020
close all; clear all
%% user defined folder inputs
figuresFolder='F:\GRABS_Data\Analyzed_SVDMethodPatch14\Figures\Hemodynamics\mCherryStateSpontCorrelationsLowFiltered';
inputFolder='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\mCherry';
%MiceAnalyze=[{'grabAM05\imaging with 575 excitation\'},{'grabAM06\imaging with 575 excitation\'},{'grabAM07\imaging with 575 excitation\'},{'grabAM08\imaging with 575 excitation\'},{'grabAM09\imaging with 575 excitation\'},{'grabAM10\imaging with 575 excitation\'}];
%MiceAnalyze=[{'grabAM09\imaging with 575 excitation'},{'grabAM08\imaging with 575 excitation'},{'SL209'},{'SL231'},{'SL232'}]; %grabs only with pre and post -physostigmine drug
%MiceAnalyze=[{'grabAM09\imaging with 575 excitation'},{'grabAM08\imaging with 575 excitation'}];%dual grabs/RCaMP mice with pre and post- physostigmine drug
%MiceAnalyze=[{'EGFP control'},{'SL224_egfp'},{'SL225_egfp'}]; % gfp control sessions
MiceAnalyze=[{'grabAM11_mcherry'},{'SL223_mcherry'}]; % mcherry control sessions
Condition='NoDrug'; %'NoDrug','PreDrug','PostDrug' ol'}]; % whethere it's a drug or drug free session
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
params.minSitDuration=5;%minimum sit duration during quiescnece state
params.CorrLags=5; %lags for cross-correlation in frames, so approximately +/- 0.5s
params.ShuffleNum=1; %number of times to shuffle data

%% add functions
addpath(genpath('F:\Sweyta\GRABS_Data\FinalPaperCodes\ParcellCorrelation'));
if ~exist(figuresFolder),mkdir(figuresFolder); end

%% get parcell indices for parcells we care about in the left hemisphere only
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

%% for each mouse, load spont and airpuff folders and perform correlations
for animal=1:length(MiceAnalyze)
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
    
    if isempty(DirFolders),continue, end
    %for each subfolder
    for folder =1:length(DirFolders)
        currfolder=fullfile(mainDir1, DirFolders(folder).name);
        %load imaging time series
        load(fullfile(currfolder,'final_dFoF_parcels.mat'),'dFoF_parcells');
        %low pass filter
        lowPassBand=[0 1]; 
        names=fieldnames(dFoF_parcells);
        for i=1:length(names)
            %filter data 
            data_filtL.(names{i})=nan(size(dFoF_parcells.(names{i}),1),size(dFoF_parcells.(names{i}),2)); 
            data_filtH.(names{i})=nan(size(dFoF_parcells.(names{i}),1),size(dFoF_parcells.(names{i}),2)); 
            for k=1:size(dFoF_parcells.(names{i}),1)
                currParcell=dFoF_parcells.(names{i})(k,:);
                if any(~isnan(currParcell))
                data_filtL.(names{i})(k,:)=eegfilt(currParcell,params.fsimaging,lowPassBand(1),lowPassBand(2));
                data_filtH.(names{i})(k,:)=currParcell-data_filtL.(names{i})(k,:);
                end 
            end
        end
        dFoF_parcells=data_filtL; 
        
        % load spike2 and pupil/face data
        load(fullfile(currfolder,'smrx_signals.mat'),'channels_data','timestamps','timing')
        files=dir(currfolder);
        for t=1:length(files)
            if contains(files(t).name,'proc','IgnoreCase',true)
                files=files(t);break;
            end
        end
%                 load(fullfile(currfolder,files.name),'proc')
%                 pupil_Norm=proc.output.pupilNorm;
%                 face_Norm=proc.output.facePC1CorrNorm;
        wheel_speed = channels_data.wheelspeed;
        
        % get pupil wheel and imaging times
        wheel_time = (1:length(wheel_speed))/params.fsspike2;
        imaging_time = timestamps.timaging;
%                 pupil_time = timing.pupilcamstart(1:length(pupil_Norm));
        
        %% extract imaging in selected parcells
        numShuffles=params.ShuffleNum;% number of random surrogate control data to build
        names=fieldnames(dFoF_parcells);
        for i=1:length(names)
            dFoF_parcells.(names{i})=dFoF_parcells.(names{i}) (CombinedParcellIdx,:); % extract data from specific parcells only
            %for each parcell, create shuffled control data by building 1000 surrogate datasets. To do this, data are randomly partionted into two slices and the order of slices are rearranged.
            dFoF_parcells_Shuffled.(names{i})=nan(size(dFoF_parcells.(names{i}),1),size(dFoF_parcells.(names{i}),2),numShuffles);
            for par=1:size(dFoF_parcells.(names{i}),1)
                for numshuff=1:numShuffles
                    data=dFoF_parcells.(names{i})(par,:);
                    frames=10:length(data);
                    splitpoint=frames(randi(length(frames)));
                    tmp1=data(1,1:splitpoint);  tmp2=data(1,splitpoint+1:end);
                    final_data=cat(2,tmp2,tmp1); %break the data into two halves and rearrange the order
                    dFoF_parcells_Shuffled.(names{i})(par,:,numshuff)=final_data;
                end
            end
        end
        %% make output figures folder
        indivFigureFolder=fullfile(figuresFolder,MiceAnalyze{animal},DirFolders(folder).name);
        if ~exist(indivFigureFolder,'dir'),mkdir(indivFigureFolder), end
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
        
%                 %% do change point detection on pupil and face to get pupil high/low arousal or face high/low movement times during sustained quiescence state
%                 %get Z-thresholds based on pupil data during quiescence, when mouse
%                 %isn't moving and when aripuffs are not given
%                 b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
%                 pupilTime_Idx=cell(1,length(sitOn_int));
%                 for st=1:length(sitOn_int)
%                     pupilTime_Idx{st}=find(pupil_time>sitOn_int(st) & pupil_time <sitOff_int(st));
%                 end
%                 pupilTime_quiescence=cell2mat(pupilTime_Idx');
%                 pupil_qui_times=pupil_time(pupilTime_quiescence);
%                 pupil_quiescence=pupil_Norm(pupilTime_quiescence);
%                 zthres_High=quantile(pupil_quiescence,0.60);
%                 zthres_Low=quantile(pupil_quiescence,0.40);
%         
%                 %get on and off timestamps for high and low arousal based on pupil
%                 [h1,Pupil_HighArousal_OnTStamp,Pupil_HighArousal_OffTStamp ] =changepoints(pupil_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
%                 title('PupilHighArousal'); saveas(h1,fullfile(indivFigureFolder,'PupilHighArousal'));
%                 [h2,Pupil_LowArousal_OnTStamp,Pupil_LowArousal_OffTStamp ] =changepoints(-pupil_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
%                 title('PupilLowArousal');saveas(h2,fullfile(indivFigureFolder,'PupilLowArousal'));
%                 %get z thresholds based on face data during quiescence only
%                 b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
%                 face_quiescence=face_Norm(pupilTime_quiescence);
%                 zthres_High=quantile(face_quiescence,0.60);
%                 zthres_Low=quantile(face_quiescence,0.40);
%         
%                 %get on and off timestamps for high and low face movment
%                 [h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
%                 title('FaceHighArousal');saveas(h3,fullfile(indivFigureFolder,'FaceHighArousal'));
%                 [h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
%                 title('FaceLowArousal');saveas(h4,fullfile(indivFigureFolder,'FaceLowArousal'));
%         
%                 % determine that pupil and face high/low arousal/movement time are at least minimum criterion seconds long
%                 idx1=find((Pupil_HighArousal_OffTStamp-Pupil_HighArousal_OnTStamp)>=params.minArousalDuration);
%                 Pupil_HighArousal_OnT_int=Pupil_HighArousal_OnTStamp(idx1); Pupil_HighArousal_OffT_int=Pupil_HighArousal_OffTStamp(idx1);
%                 idx2=find((Pupil_LowArousal_OffTStamp-Pupil_LowArousal_OnTStamp)>=params.minArousalDuration);
%                 Pupil_LowArousal_OnT_int=Pupil_LowArousal_OnTStamp(idx2); Pupil_LowArousal_OffT_int=Pupil_LowArousal_OffTStamp(idx2);
%                 idx3=find((Face_HighArousal_OffTStamp-Face_HighArousal_OnTStamp)>=params.minArousalDuration);
%                 Face_HighArousal_OnT_int=Face_HighArousal_OnTStamp(idx3); Face_HighArousal_OffT_int=Face_HighArousal_OffTStamp(idx3);
%                 idx4=find((Face_LowArousal_OffTStamp-Face_LowArousal_OnTStamp)>=params.minArousalDuration);
%                 Face_LowArousal_OnT_int=Face_LowArousal_OnTStamp(idx4); Face_LowArousal_OffT_int=Face_LowArousal_OffTStamp(idx4);
%         
%                 % get pupil on and face on times if both on and off times occur entirely during sustained queiscence states identified in the previous step
%                 toDelete=ones(1,length(Pupil_HighArousal_OnT_int));
%                 for rj=1:length(Pupil_HighArousal_OnT_int)
%                     tmp = find (Pupil_HighArousal_OnT_int(rj)>=sitOn_final & Pupil_HighArousal_OffT_int(rj)<=sitOff_final);
%                     toDelete(rj)=isempty(tmp);
%                 end
%                 Pupil_HighArousal_On_final=Pupil_HighArousal_OnT_int(~toDelete);
%                 Pupil_HighArousal_Off_final=Pupil_HighArousal_OffT_int(~toDelete);
%         
%                 toDelete=ones(1,length(Pupil_LowArousal_OnT_int));
%                 for rj=1:length(Pupil_LowArousal_OnT_int)
%                     tmp = find (Pupil_LowArousal_OnT_int(rj)>=sitOn_final & Pupil_LowArousal_OffT_int(rj)<=sitOff_final);
%                     toDelete(rj)=isempty(tmp);
%                 end
%                 Pupil_LowArousal_On_final=Pupil_LowArousal_OnT_int(~toDelete);
%                 Pupil_LowArousal_Off_final=Pupil_LowArousal_OffT_int(~toDelete);
%         
%                 toDelete=ones(1,length(Face_HighArousal_OnT_int));
%                 for rj=1:length(Face_HighArousal_OnT_int)
%                     tmp = find (Face_HighArousal_OnT_int(rj)>=sitOn_final & Face_HighArousal_OffT_int(rj)<=sitOff_final);
%                     toDelete(rj)=isempty(tmp);
%                 end
%                 Face_HighArousal_On_final=Face_HighArousal_OnT_int(~toDelete);
%                 Face_HighArousal_Off_final=Face_HighArousal_OffT_int(~toDelete);
%         
%                 toDelete=ones(1,length(Face_LowArousal_OnT_int));
%                 for rj=1:length(Face_LowArousal_OnT_int)
%                     tmp = find (Face_LowArousal_OnT_int(rj)>=sitOn_final & Face_LowArousal_OffT_int(rj)<=sitOff_final);
%                     toDelete(rj)=isempty(tmp);
%                 end
%                 Face_LowArousal_On_final=Face_LowArousal_OnT_int(~toDelete);
%                 Face_LowArousal_Off_final=Face_LowArousal_OffT_int(~toDelete);
%         
        %% plot on and off times for sanity check
        %plot locomotion state on off times
        h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified states
        set(0,'CurrentFigure',h);ax1=subplot(6,1,1);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('LocomotionState');
        for tt=1:length(wheelOn_final)
            plot([wheelOn_final(tt),wheelOn_final(tt)], ylimits, 'g');
            plot([wheelOff_final(tt),wheelOff_final(tt)], ylimits, 'r');
        end
        %plot airpuff timea andn imaging onset/offset
        if ~isempty(allEvts)
            for tt=1:length(allEvts)
                plot([allEvts(tt),allEvts(tt)], ylimits, 'k');
            end
        end
        plot([imaging_time(1),imaging_time(1)], ylimits, 'm');
        plot([imaging_time(end),imaging_time(end)], ylimits, 'm');
        
        %plot quiescence state on off times
        set(0,'CurrentFigure',h);ax2=subplot(6,1,2);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('QuiescenceState');
        for tt=1:length(sitOn_final)
            plot([sitOn_final(tt),sitOn_final(tt)], ylimits, 'g');
            plot([sitOff_final(tt),sitOff_final(tt)], ylimits, 'r');
        end
        %plot airpuff timea and imaging onset/offset times
        if exist('allEvts','var')
            for tt=1:length(allEvts)
                plot([allEvts(tt),allEvts(tt)], ylimits, 'k');
            end
        end
        plot([imaging_time(1),imaging_time(1)], ylimits, 'm');
        plot([imaging_time(end),imaging_time(end)], ylimits, 'm');
        
%                 %plot pupil and face state on off times
%                 set(0,'CurrentFigure',h);ax3=subplot(6,1,3);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilHighState');
%                 for tt=1:length(Pupil_HighArousal_On_final)
%                     plot([Pupil_HighArousal_On_final(tt),Pupil_HighArousal_On_final(tt)], ylimits, 'g');
%                     plot([Pupil_HighArousal_Off_final(tt),Pupil_HighArousal_Off_final(tt)], ylimits, 'r');
%                 end
%         
%                 set(0,'CurrentFigure',h);ax4=subplot(6,1,4);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilLowState');
%                 for tt=1:length(Pupil_LowArousal_On_final)
%                     plot([Pupil_LowArousal_On_final(tt),Pupil_LowArousal_On_final(tt)], ylimits, 'g');
%                     plot([Pupil_LowArousal_Off_final(tt),Pupil_LowArousal_Off_final(tt)], ylimits, 'r');
%                 end
%         
%                 set(0,'CurrentFigure',h);ax5=subplot(6,1,5);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceHighState');
%                 for tt=1:length(Face_HighArousal_On_final)
%                     plot([Face_HighArousal_On_final(tt),Face_HighArousal_On_final(tt)], ylimits, 'g');
%                     plot([Face_HighArousal_Off_final(tt),Face_HighArousal_Off_final(tt)], ylimits, 'r');
%                 end
%         
%                 set(0,'CurrentFigure',h);ax6=subplot(6,1,6);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceLowState');
%                 for tt=1:length(Face_LowArousal_On_final)
%                     plot([Face_LowArousal_On_final(tt),Face_LowArousal_On_final(tt)], ylimits, 'g');
%                     plot([Face_LowArousal_Off_final(tt),Face_LowArousal_Off_final(tt)], ylimits, 'r');
%                 end
%                 linkaxes([ax1, ax2, ax3,ax4,ax5,ax6],'x');
        saveas(h,fullfile(indivFigureFolder,'BehavioralStatesFinalOnOffTimes'));
        %% ensure for comparison purpose, wihitn each session, the number of sit and run trials as well as the total time within each trial is matched. If no wheel or sit trials exist, skip
        if ~isempty(wheelOn_final) && ~isempty(sitOn_final)
            numWheel=length(wheelOn_final); numSit=length(sitOn_final);
            totalTrials_wheel=min(numWheel,numSit);
            wheelOn_final1=wheelOn_final(1:totalTrials_wheel);
            wheelOff_final1=wheelOff_final(1:totalTrials_wheel);
            sitOn_final1=sitOn_final(1:totalTrials_wheel);
            sitOff_final1=sitOff_final(1:totalTrials_wheel);
            
            wheelOnDur=wheelOff_final1-wheelOn_final1;
            sitOnDur=sitOff_final1-sitOn_final1;
            minDur=min(wheelOnDur,sitOnDur);
            wheelOff_final1=wheelOn_final1+minDur;
            sitOff_final1=sitOn_final1+minDur;
            
            % run cross correlations on sit and run data and extract the max correlation
            [loccorrMat{animal,folder},loccorrMat_Shuffled{animal,folder},sitcorrMat{animal,folder},sitcorrMat_Shuffled{animal,folder},locImaging{animal,folder},sitImaging{animal,folder},locImaging_Shuffled{animal,folder},sitImaging_Shuffled{animal,folder}]...
                =stateCrossCorr(wheelOn_final1,wheelOff_final1,sitOn_final1,sitOff_final1,imaging_time, dFoF_parcells,dFoF_parcells_Shuffled,params.CorrLags,names,params.signalsExtraction);
        end
%         %% ensure for comparison purpose, wihitn each session, the high and low pupil arousal trials as well as the total time within each trial is matched. If no high/low arousal trials exist, skip
%                 if ~isempty(Pupil_HighArousal_On_final) && ~isempty(Pupil_LowArousal_On_final)
%                     numHPupil=length(Pupil_HighArousal_On_final); numLPupil=length(Pupil_LowArousal_On_final);
%                     totalTrials_Pupil=min(numHPupil,numLPupil);
%                     Pupil_HighArousal_On_final1=Pupil_HighArousal_On_final(1:totalTrials_Pupil);
%                     Pupil_HighArousal_Off_final1=Pupil_HighArousal_Off_final(1:totalTrials_Pupil);
%                     Pupil_LowArousal_On_final1=Pupil_LowArousal_On_final(1:totalTrials_Pupil);
%                     Pupil_LowArousal_Off_final1=Pupil_LowArousal_Off_final(1:totalTrials_Pupil);
%         
%                     Pupil_HighArousal_OnDur=Pupil_HighArousal_Off_final1-Pupil_HighArousal_On_final1;
%                     Pupil_LowArousal_OnDur=Pupil_LowArousal_Off_final1-Pupil_LowArousal_On_final1;
%                     minDur=min(Pupil_HighArousal_OnDur,Pupil_LowArousal_OnDur);
%                     Pupil_HighArousal_Off_final1=Pupil_HighArousal_On_final1+minDur;
%                     Pupil_LowArousal_Off_final1=Pupil_LowArousal_On_final1+minDur;
%         
%                     %run cross corelations on high and low pupil trials
%                     [PupilHighcorrMat{animal,folder},PupilHighcorrMat_Shuffled{animal,folder},PupilLowcorrMat{animal,folder},PupilLowcorrMat_Shuffled{animal,folder},PupilHighImaging{animal,folder},PupilLowImaging{animal,folder},PupilHighImaging_Shuffled{animal,folder},PupilLowImaging_Shuffled{animal,folder}]...
%                         =stateCrossCorr(Pupil_HighArousal_On_final1,Pupil_HighArousal_Off_final1,Pupil_LowArousal_On_final1,Pupil_LowArousal_Off_final1,imaging_time, dFoF_parcells,dFoF_parcells_Shuffled,params.CorrLags,names,params.signalsExtraction);
%                 end
%         
%                 %% ensure for comparison purpose, wihitn each session, the number of high and low facial motion trials as well as the total time within each trial is matched. If no high/low facial motion trials exist, skip
%                 if ~isempty(Face_HighArousal_On_final) && ~isempty(Face_LowArousal_On_final)
%                     numFaceH=length(Face_HighArousal_On_final); numFaceL=length(Face_LowArousal_On_final);
%                     totalTrials_Face=min(numFaceH,numFaceL);
%                     Face_HighArousal_On_final1=Face_HighArousal_On_final(1:totalTrials_Face);
%                     Face_HighArousal_Off_final1=Face_HighArousal_Off_final(1:totalTrials_Face);
%                     Face_LowArousal_On_final1=Face_LowArousal_On_final(1:totalTrials_Face);
%                     Face_LowArousal_Off_final1=Face_LowArousal_Off_final(1:totalTrials_Face);
%         
%                     Face_HighArousal_OnDur=Face_HighArousal_Off_final1-Face_HighArousal_On_final1;
%                     Face_LowArousal_OnDur=Face_LowArousal_Off_final1-Face_LowArousal_On_final1;
%                     minDur=min(Face_HighArousal_OnDur,Face_LowArousal_OnDur);
%                     Face_HighArousal_Off_final1=Face_HighArousal_On_final1+minDur;
%                     Face_LowArousal_Off_final1=Face_LowArousal_On_final1+minDur;
%         
%                     %run correlations on non-filtered, low and high filtered data
%                     [FaceHighcorrMat{animal,folder},FaceHighcorrMat_Shuffled{animal,folder},FaceLowcorrMat{animal,folder},FaceLowcorrMat_Shuffled{animal,folder},FaceHighImaging{animal,folder},FaceLowImaging{animal,folder},FaceHighImaging_Shuffled{animal,folder},FaceLowImaging_Shuffled{animal,folder}]...
%                         =stateCrossCorr(Face_HighArousal_On_final1,Face_HighArousal_Off_final1,Face_LowArousal_On_final1,Face_LowArousal_Off_final1,imaging_time, dFoF_parcells,dFoF_parcells_Shuffled,params.CorrLags,names,params.signalsExtraction);
%                 end
    end
    %average correlation matrices across trials concatenated across sessions within an animal
    indivFigureFolder=fullfile(figuresFolder,MiceAnalyze{animal});
    names=fieldnames(dFoF_parcells);
    if strcmp(params.signalsExtraction,'RCaMP_AC'),names=[names;{'BG'}]; end
    for i=1:length(names)
        if  size(loccorrMat,1)==animal
            [hs,loccorrMatCatZScore.(names{i}){animal},sitcorrMatCatZScore.(names{i}){animal},loccorrMatCatZScoreShuffled.(names{i}){animal},sitcorrMatCatZScoreShuffled.(names{i}){animal}...
                ,loccorrMatCat.(names{i}){animal},sitcorrMatCat.(names{i}){animal},loccorrMatCatShuffled.(names{i}){animal}, sitcorrMatCatShuffled.(names{i}){animal}]...
                = averageCorrMatrixSessions(loccorrMat(animal,:),sitcorrMat(animal,:),loccorrMat_Shuffled(animal,:),sitcorrMat_Shuffled(animal,:),names{i},'Locomotion','Quiescence',parcellnames);
            
            saveas(hs,fullfile(indivFigureFolder,strcat(names{i},'LocomotionvsQuiescenceCorrMatrix')));
        end
%         
%                 if size(PupilHighcorrMat,1)==animal
%                     [hs,PupilHighcorrMatCatZScore.(names{i}){animal},PupilLowcorrMatCatZScore.(names{i}){animal},PupilHighcorrMatCatZScoreShuffled.(names{i}){animal},PupilLowcorrMatCatZScoreShuffled.(names{i}){animal}...
%                         ,PupilHighcorrMatCat.(names{i}){animal},PupilLowcorrMatCat.(names{i}){animal},PupilHighcorrMatCatShuffled.(names{i}){animal}, PupilLowcorrMatCatShuffled.(names{i}){animal}]...
%                         = averageCorrMatrixSessions(PupilHighcorrMat(animal,:),PupilLowcorrMat(animal,:),PupilHighcorrMat_Shuffled(animal,:),PupilLowcorrMat_Shuffled(animal,:),names{i},'PupilHigh','PupilLow',parcellnames);
%         
%                     saveas(hs,fullfile(indivFigureFolder,strcat(names{i},'HighvsLowPupilCorrMatrix')));
%                 end
%                 if size(FaceHighcorrMat,1)==animal
%                     [hs,FaceHighcorrMatCatZScore.(names{i}){animal},FaceLowcorrMatCatZScore.(names{i}){animal},FaceHighcorrMatCatZScoreShuffled.(names{i}){animal},FaceLowcorrMatCatZScoreShuffled.(names{i}){animal}...
%                         ,FaceHighcorrMatCat.(names{i}){animal},FaceLowcorrMatCat.(names{i}){animal},FaceHighcorrMatCatShuffled.(names{i}){animal}, FaceLowcorrMatCatShuffled.(names{i}){animal}]...
%                         = averageCorrMatrixSessions(FaceHighcorrMat(animal,:),FaceLowcorrMat(animal,:),FaceHighcorrMat_Shuffled(animal,:),FaceLowcorrMat_Shuffled(animal,:),names{i},'FaceHigh','FaceLow',parcellnames);
%         
%                     saveas(hs,fullfile(indivFigureFolder,strcat(names{i},'HighvsLowFaceCorrMatrix')));
%                 end
    end
    close all;
end
%save output
% save(fullfile(figuresFolder,'IndivMouseOutput.mat'),'loccorrMatCatZScore','sitcorrMatCatZScore','loccorrMatCatZScoreShuffled','sitcorrMatCatZScoreShuffled','PupilHighcorrMatCatZScore','PupilLowcorrMatCatZScore',...
%     'PupilHighcorrMatCatZScoreShuffled','PupilLowcorrMatCatZScoreShuffled','FaceHighcorrMatCatZScore','FaceLowcorrMatCatZScore','FaceHighcorrMatCatZScoreShuffled','FaceLowcorrMatCatZScoreShuffled');
%% population averages, take a mean across animals, first mean of fisher's z , then backtransform z to pearson's r for display purposes
names=fieldnames(loccorrMatCatZScore);
for i=1:length(names)
    %% do stats on locomotion vs quiescence
    figFolder=fullfile(figuresFolder,'Locomotion'); if ~exist(figFolder),mkdir(figFolder), end
    [output.(names{i}).Loc.meanState1,output.(names{i}).Loc.meanState2,output.(names{i}).Loc.meanDiffMat,output.(names{i}).Loc.V1S1M2]=statsCorrelation1(loccorrMatCatZScore.(names{i}),sitcorrMatCatZScore.(names{i}), names{i}, pixelpropSizeParcell,'loc','sit',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
    figFolder=fullfile(figuresFolder,'LocomotionShuffled'); if ~exist(figFolder),mkdir(figFolder), end
    [output.(names{i}).LocShuffled.meanState1,output.(names{i}).LocShuffled.meanState2,output.(names{i}).LocShuffled.meanDiffMat,output.(names{i}).LocShuffled.V1S1M2]=statsCorrelation1(loccorrMatCatZScoreShuffled.(names{i}),sitcorrMatCatZScoreShuffled.(names{i}), names{i}, pixelpropSizeParcell,'locShuff','sitShuff',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
%         %% do stats on pupil high vs low
%         figFolder=fullfile(figuresFolder,'Pupil'); if ~exist(figFolder),mkdir(figFolder), end
%         [output.(names{i}).Pupil.meanState1,output.(names{i}).Pupil.meanState2,output.(names{i}).Pupil.meanDiffMat,output.(names{i}).Pupil.V1S1M2]=statsCorrelation1(PupilHighcorrMatCatZScore.(names{i}),PupilLowcorrMatCatZScore.(names{i}), names{i}, pixelpropSizeParcell,'PupilHigh','PupilLow',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
%         figFolder=fullfile(figuresFolder,'PupilShuffled'); if ~exist(figFolder),mkdir(figFolder), end;
%         [output.(names{i}).PupilShuffled.meanState1,output.(names{i}).PupilShuffled.meanState2,output.(names{i}).PupilShuffled.meanDiffMat,output.(names{i}).PupilShuffled.V1S1M2]=statsCorrelation1(PupilHighcorrMatCatZScoreShuffled.(names{i}),PupilLowcorrMatCatZScoreShuffled.(names{i}), names{i}, pixelpropSizeParcell,'PupilHighShuff','PupilLowShuff',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
%     
%         %% do stats on face high vs low
%         figFolder=fullfile(figuresFolder,'Face'); if ~exist(figFolder),mkdir(figFolder), end
%         [output.(names{i}).Face.meanState1,output.(names{i}).Face.meanState2,output.(names{i}).Face.meanDiffMat,output.(names{i}).Face.V1S1M2]=statsCorrelation1(FaceHighcorrMatCatZScore.(names{i}),FaceLowcorrMatCatZScore.(names{i}), names{i}, pixelpropSizeParcell,'FaceHigh','FaceLow',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
%         figFolder=fullfile(figuresFolder,'FaceShuffled'); if ~exist(figFolder),mkdir(figFolder), end
%         [output.(names{i}).FaceShuffled.meanState1,output.(names{i}).FaceShuffled.meanState2,output.(names{i}).FaceShuffled.meanDiffMat,output.(names{i}).FaceShuffled.V1S1M2]=statsCorrelation1(FaceHighcorrMatCatZScoreShuffled.(names{i}),FaceLowcorrMatCatZScoreShuffled.(names{i}), names{i}, pixelpropSizeParcell,'FaceHighShuff','FaceLowShuff',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
end

%save output
save(fullfile(figuresFolder,'SummaryData.mat'),'output');


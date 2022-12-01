%% Figure 3 This script calculates EEG power during three sustained states (locomotion, face high, face low). All time are in seconds
%written by Sweyta Lohani 2020
close all; clear all
%% user defined folder inputs
figuresFolder='F:\Figures';%where figure will be saved 
inputFolder='W:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice';%where input data are located 
MiceAnalyze=[{'grabAM05\imaging with 575 excitation\'},{'grabAM06\imaging with 575 excitation\'},{'grabAM07\imaging with 575 excitation\'},{'grabAM08\imaging with 575 excitation\'},{'grabAM09\imaging with 575 excitation\'},{'grabAM10\imaging with 575 excitation\'}];
Condition='NoDrug'; %'NoDrug','PreDrug','PostDrug' ol'}]; % whethere it's a drug or drug free session
%% idx of sessions that don't have artifacts in EEG for the  sessions belonging to mice specified above 
goodEEGidx=[1 1 1 0; 0 1 1 0;1 1 1 0;1 0 0 0; 1 1 0 0; 0 0 1 0]; %manually curated sessions that don't have artifacts 
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
%% add functions
addpath(genpath('./Functions'));
if ~exist(figuresFolder),mkdir(figuresFolder); end
%% for each mouse, load spont and airpuff folders and perform correlations
for animal=1:length(MiceAnalyze)
    mainDir1 =fullfile(inputFolder,MiceAnalyze{animal});
    folders=dir(mainDir1);
    dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
    DirFolders= folders(dirFlags);
    noDrugFlag=~contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim and drug sessions
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
        if goodEEGidx(animal,folder)==0, continue, end; % skip if the EEG is bad for this session 
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
        wheel_speed = channels_data.wheelspeed(1:end-1);%remove last element which is nan
        EEG=channels_data.EEG(1:end-1);%raw EEG, remove last element which is nan
        
        %filter EEG signal at 1-10 Hz and get amplitude envelope 
        [filteredEEG] = eegfilt(EEG',params.fsspike2,1,10); 
        EEG_env=abs(hilbert(filteredEEG));
        
        % get pupil wheel and imaging times
        wheel_time = (1:length(wheel_speed))/params.fsspike2;
        EEG_time=(1:length(EEG))/params.fsspike2;
        imaging_time = timestamps.timaging;
        pupil_time = timing.pupilcamstart(1:length(pupil_Norm));
        
        %% get locomotion on/off and quiescence on off times
        %locomotion periods should be at least some criterion s long with some criterion sec since locomotion onset, some criterion sec before locomotion offset, excluding any events (airpuff/stim)
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
        
        %% do change point detection on face to get face high/low movement times during sustained quiescence state
        %get Z-thresholds based on face data during quiescence, when mouse isn't moving and when aripuffs are not given
        b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
        pupilTime_Idx=cell(1,length(sitOn_int));
        for st=1:length(sitOn_int)
            pupilTime_Idx{st}=find(pupil_time>sitOn_int(st) & pupil_time <sitOff_int(st));
        end
        pupilTime_quiescence=cell2mat(pupilTime_Idx');
        face_quiescence=face_Norm(pupilTime_quiescence);
        zthres_High=quantile(face_quiescence,0.60);
        zthres_Low=quantile(face_quiescence,0.40);
        
        %get on and off timestamps for high and low face movment
        [h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
        title('FaceHighArousal');
        [h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
        title('FaceLowArousal');
        
        % determine that face high/low arousal/movement time are at least minimum criterion seconds long
        idx3=find((Face_HighArousal_OffTStamp-Face_HighArousal_OnTStamp)>=params.minArousalDuration);
        Face_HighArousal_OnT_int=Face_HighArousal_OnTStamp(idx3); Face_HighArousal_OffT_int=Face_HighArousal_OffTStamp(idx3);
        idx4=find((Face_LowArousal_OffTStamp-Face_LowArousal_OnTStamp)>=params.minArousalDuration);
        Face_LowArousal_OnT_int=Face_LowArousal_OnTStamp(idx4); Face_LowArousal_OffT_int=Face_LowArousal_OffTStamp(idx4);
        
        % get  face on times if both on and off times occur entirely during sustained queiscence states identified in the previous step       
        toDelete=ones(1,length(Face_HighArousal_OnT_int));
        for rj=1:length(Face_HighArousal_OnT_int)
            tmp = find (Face_HighArousal_OnT_int(rj)>=sitOn_final & Face_HighArousal_OffT_int(rj)<=sitOff_final);
            toDelete(rj)=isempty(tmp);
        end
        Face_HighArousal_On_final=Face_HighArousal_OnT_int(~toDelete);
        Face_HighArousal_Off_final=Face_HighArousal_OffT_int(~toDelete);
        
        toDelete=ones(1,length(Face_LowArousal_OnT_int));
        for rj=1:length(Face_LowArousal_OnT_int)
            tmp = find (Face_LowArousal_OnT_int(rj)>=sitOn_final & Face_LowArousal_OffT_int(rj)<=sitOff_final);
            toDelete(rj)=isempty(tmp);
        end
        Face_LowArousal_On_final=Face_LowArousal_OnT_int(~toDelete);
        Face_LowArousal_Off_final=Face_LowArousal_OffT_int(~toDelete);
        
        
      %% ensure for comparison purpose, wihitn each session, the number of face low (sit), face high (sit), and locomotion trials are matched. If no wheel or face high or face low trials exist, skip
        if ~isempty(wheelOn_final) &&~isempty(Face_HighArousal_On_final) && ~isempty(Face_LowArousal_On_final)
            %sort epochs by duration and use the same number of epochs across states 
            wheelOnDur=(wheelOff_final-wheelOn_final)';
            Face_HighDur=Face_HighArousal_Off_final-Face_HighArousal_On_final;
            Face_LowDur=Face_LowArousal_Off_final-Face_LowArousal_On_final;
            numWheel=length(wheelOn_final); numFaceH=length(Face_HighArousal_On_final); numFaceL=length(Face_LowArousal_On_final);
            totalTrials_wheel=min([numWheel,numFaceH,numFaceL]);
            [sortedWheelDur,sortedWheelidx]=sort(wheelOnDur,'descend');  
            [sortedFaceHighDur,sortedFaceHighidx]=sort(Face_HighDur,'descend');
            [sortedFaceLowDur,sortedFaceLowidx]=sort(Face_LowDur,'descend');
            sortedWheelIdx=sortedWheelidx(1:totalTrials_wheel); 
            sortedFaceHighIdx=sortedFaceHighidx(1:totalTrials_wheel); 
            sortedFaceLowIdx=sortedFaceLowidx(1:totalTrials_wheel); 
            
            wheelOn_final1=wheelOn_final(sortedWheelIdx);
            wheelOff_final1=wheelOff_final(sortedWheelIdx);
            Face_HighArousal_On_final1=Face_HighArousal_On_final(sortedFaceHighIdx);
            Face_HighArousal_Off_final1=Face_HighArousal_Off_final(sortedFaceHighIdx);
            Face_LowArousal_On_final1=Face_LowArousal_On_final(sortedFaceLowIdx);
            Face_LowArousal_Off_final1=Face_LowArousal_Off_final(sortedFaceLowIdx);
            
            %match the duration of epochs across states 
            wheelOnDur=(wheelOff_final1-wheelOn_final1)';
            Face_HighArousal_OnDur=Face_HighArousal_Off_final1-Face_HighArousal_On_final1;
            Face_LowArousal_OnDur=Face_LowArousal_Off_final1-Face_LowArousal_On_final1;
            minDur=min([Face_HighArousal_OnDur,Face_LowArousal_OnDur,wheelOnDur],[],2);
            wheelOff_final1=wheelOn_final1'+minDur;
            Face_HighArousal_Off_final1=Face_HighArousal_On_final1+minDur;
            Face_LowArousal_Off_final1=Face_LowArousal_On_final1+minDur;
            
            CompData.loc.TotalDur{animal,folder}=sum(minDur);
            CompData.loc.numEpochs{animal,folder}=totalTrials_wheel;
            
            %extract EEG and bandpower during states 
            [loc.EEG{animal,folder},faceHigh.EEG{animal,folder},faceLow.EEG{animal,folder},loc.EEG_env{animal,folder},faceHigh.EEG_env{animal,folder},faceLow.EEG_env{animal,folder},...
                loc.lowpower{animal,folder},faceHigh.lowpower{animal,folder},faceLow.lowpower{animal,folder},loc.highpower{animal,folder},faceHigh.highpower{animal,folder},faceLow.highpower{animal,folder}]...
                =stateEEGAnalyzeCombined(wheelOn_final1,wheelOff_final1,Face_HighArousal_On_final1,Face_HighArousal_Off_final1,Face_LowArousal_On_final1,Face_LowArousal_Off_final1,EEG_time, EEG,EEG_env,params.fsspike2);
 
        end       
    end
 % concatenate data across sessions and average to get one value per animal 
 if sum(goodEEGidx(animal,:))>0 %if there is at least one good EEG session for the animal, otherwise skip 
    loc.highpower_animal{animal}=nanmean(cell2mat(loc.highpower(animal,:))); 
    loc.lowpower_animal{animal}=nanmean(cell2mat(loc.lowpower(animal,:))); 
    loc.highLowRatio_animal{animal}=nanmean(cell2mat(loc.highpower(animal,:))./cell2mat(loc.lowpower(animal,:))); 
    loc.lowAmpEnv_animal{animal}=nanmean(cell2mat(loc.EEG_env(animal,:))); 

    FaceH.highpower_animal{animal}=nanmean(cell2mat(faceHigh.highpower(animal,:))); 
    FaceL.highpower_animal{animal}=nanmean(cell2mat(faceLow.highpower(animal,:))); 
    FaceH.lowpower_animal{animal}=nanmean(cell2mat(faceHigh.lowpower(animal,:))); 
    FaceL.lowpower_animal{animal}=nanmean(cell2mat(faceLow.lowpower(animal,:))); 
    FaceH.highLowRatio_animal{animal}=nanmean(cell2mat(faceHigh.highpower(animal,:))./cell2mat(faceHigh.lowpower(animal,:))); 
    FaceL.highLowRatio_animal{animal}=nanmean(cell2mat(faceHigh.highpower(animal,:))./cell2mat(faceLow.lowpower(animal,:)));
    FaceH.lowAmpEnv_animal{animal}=nanmean(cell2mat(faceHigh.EEG_env(animal,:))); 
    FaceL.lowAmpEnv_animal{animal}=nanmean(cell2mat(faceLow.EEG_env(animal,:)));  
 end    
end

%% make plots of high and low power and do statistics 
figure1=figure;
subplot(2,1,1); title('HighPower'); hold on; 
data=[cell2mat(FaceL.highpower_animal)',cell2mat(FaceH.highpower_animal)',cell2mat(loc.highpower_animal)']; 
stats.highpower.data=data; 
%do one way repeated measures anova followed by posthoc t-tests
t = table(data(:,1),data(:,2),data(:,3),...
            'VariableNames',{'State1','State2','State3'});
Method = table([1 2 3]','VariableNames',{'States'});
RM = fitrm(t,'State1-State3~1','WithinDesign',Method);
stats.highpower.RM= ranova(RM);
[stats.highpower.posthoc.face_h, stats.highpower.posthoc.face_p,~,stats.highpower.posthoc.face_stats]=ttest(data(:,1),data(:,2)); %paired ttest between low and high face
[stats.highpower.posthoc.loc_h, stats.highpower.posthoc.loc_p,~,stats.highpower.posthoc.loc_stats]=ttest(data(:,2),data(:,3)); %paired tttest between high face and locomotion 
plot(data','-ko','MarkerFaceColor','k'); hold on;
currMean=mean(data,1);
plot(currMean,'+r','MarkerSize',30);  %draw markers at mean
xlim([0.5 3.5]);xticks(1:3); xticklabels([{'faceL'},{'faceH'},{'loc'}]); xtickangle(90); box off; ylabel('HighPower'); ylim([0 0.0015]); set(gca,'TickDir','out'); 

subplot(2,1,2); title('LowPower'); hold on; 
data=[cell2mat(FaceL.lowpower_animal)',cell2mat(FaceH.lowpower_animal)',cell2mat(loc.lowpower_animal)']; 
stats.lowpower.data=data; 
%do one way repeated measures anova followed by posthoc t-tests
t = table(data(:,1),data(:,2),data(:,3),...
            'VariableNames',{'State1','State2','State3'});
Method = table([1 2 3]','VariableNames',{'States'});
RM = fitrm(t,'State1-State3~1','WithinDesign',Method);
stats.lowpower.RM= ranova(RM);
[stats.lowpower.posthoc.face_h, stats.lowpower.posthoc.face_p,~,stats.lowpower.posthoc.face_stats]=ttest(data(:,1),data(:,2)); %paired ttest between low and high face
[stats.lowpower.posthoc.loc_h, stats.lowpower.posthoc.loc_p,~,stats.lowpower.posthoc.loc_stats]=ttest(data(:,2),data(:,3)); %paired tttest between high face and locomotion 
plot(data','-ko','MarkerFaceColor','k'); hold on;
currMean=mean(data,1);
plot(currMean,'+r','MarkerSize',30);  %draw markers at mean
xlim([0.5 3.5]);xticks(1:3); xticklabels([{'faceL'},{'faceH'},{'loc'}]); xtickangle(90); box off; ylabel('LowPower'); ylim([0 0.015]); set(gca,'TickDir','out'); 

save(fullfile(figuresFolder,'Summary.mat'),'stats'); 
saveas(figure1,fullfile(figuresFolder,'EEGSustainedStates')); saveas(figure1,fullfile(figuresFolder,'EEGSustainedStates.pdf')); 



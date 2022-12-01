function[wheelOn_final,wheelOff_final,Face_HighArousal_On_final,Face_HighArousal_Off_final,Face_LowArousal_On_final,Face_LowArousal_Off_final]...
    =getSustainedStates(timing,params,imaging_time,pupil_time,wheel_time, face_Norm,wheel_speed, files,indivFigureFolder)
%this function gets sustained locomotion, face high, and face low states
%written by SL, 2022
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
sitOn=[0;wheelOff]; %use 0 as the first sit on time;
sitOff=[wheelOn;imaging_time(end)];%use wheelOn times as sit off times;


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
if length(files) ==1
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
    title('FaceHighArousal');saveas(h3,fullfile(indivFigureFolder,'FaceHighArousal'));
    [h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
    title('FaceLowArousal');saveas(h4,fullfile(indivFigureFolder,'FaceLowArousal'));
    
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
    Face_HighArousal_On_final=Face_HighArousal_OnT_int(~toDelete)';
    Face_HighArousal_Off_final=Face_HighArousal_OffT_int(~toDelete)';
    
    toDelete=ones(1,length(Face_LowArousal_OnT_int));
    for rj=1:length(Face_LowArousal_OnT_int)
        tmp = find (Face_LowArousal_OnT_int(rj)>=sitOn_final & Face_LowArousal_OffT_int(rj)<=sitOff_final);
        toDelete(rj)=isempty(tmp);
    end
    Face_LowArousal_On_final=Face_LowArousal_OnT_int(~toDelete)';
    Face_LowArousal_Off_final=Face_LowArousal_OffT_int(~toDelete)';
else
    Face_HighArousal_On_final=[];
    Face_HighArousal_Off_final=[];
    Face_LowArousal_On_final=[];
    Face_LowArousal_Off_final=[];
end
if params.plotStates
    %% plot on and off times for sanity check
    %plot locomotion state on off times
    h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified states
    set(0,'CurrentFigure',h);ax1=subplot(4,1,1);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('LocomotionState');
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
    set(0,'CurrentFigure',h);ax2=subplot(4,1,2);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('QuiescenceState');
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
    if length(files) ==1
        %plot face state on off times
        set(0,'CurrentFigure',h);ax5=subplot(4,1,3);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceHighState');
        for tt=1:length(Face_HighArousal_On_final)
            plot([Face_HighArousal_On_final(tt),Face_HighArousal_On_final(tt)], ylimits, 'g');
            plot([Face_HighArousal_Off_final(tt),Face_HighArousal_Off_final(tt)], ylimits, 'r');
        end
        
        set(0,'CurrentFigure',h);ax6=subplot(4,1,4);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceLowState');
        for tt=1:length(Face_LowArousal_On_final)
            plot([Face_LowArousal_On_final(tt),Face_LowArousal_On_final(tt)], ylimits, 'g');
            plot([Face_LowArousal_Off_final(tt),Face_LowArousal_Off_final(tt)], ylimits, 'r');
        end
    end
    linkaxes([ax1, ax2,ax5,ax6],'x');
    saveas(h,fullfile(indivFigureFolder,'BehavioralStatesFinalOnOffTimes'));
end
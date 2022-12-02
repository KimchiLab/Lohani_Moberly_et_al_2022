function [tStamps]=extractFacePupilTStamps(params,channels_data,proc,timestamps, timing,indivFigureFolder)
%% This function extract pupil and face state transition (transition from low to high arousal) on/off timestamps 
%Sweyta Lohani, 2020 
%% load state data and times 
%get face pupil and running data 
face_Norm=proc.output.facePC1CorrNorm;
%pupil_Norm=proc.output.pupilNorm;
wheel_speed = channels_data.wheelspeed;

% get pupil wheel and imaging times
wheel_time = (1:length(wheel_speed))/params.fsspike2;
imaging_time = timestamps.timaging;
pupil_time = timing.pupilcamstart(1:length(face_Norm));

% get locomotion on/off and quiescence on off times
wheelOn=timing.allwheelon;
wheelOff=timing.allwheeloff;

%% queiscence should be at least some criterion s long with some criterion s since locomotion offset
%and some criterion s before subsequent locomotion onset,excluding any events (airpuff/visual stim)
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

%% do change point detection on face PC1
% b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
pupilTime_Idx=cell(1,length(sitOn_int));
for st=1:length(sitOn_int)
    pupilTime_Idx{st}=find(pupil_time>sitOn_int(st) & pupil_time <sitOff_int(st));
end
pupilTime_quiescence=cell2mat(pupilTime_Idx');
% pupil_quiescence=pupil_Norm(pupilTime_quiescence);
% zthres_High=quantile(pupil_quiescence,0.60);
% zthres_Low=quantile(pupil_quiescence,0.40);

% %get on and off timestamps for high and low arousal based on pupil
% [h1,Pupil_HighArousal_OnTStamp,Pupil_HighArousal_OffTStamp ] =changepoints(pupil_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
% title('PupilHighArousal'); saveas(h1,fullfile(indivFigureFolder,'PupilHighArousal'));
% [h2,Pupil_LowArousal_OnTStamp,Pupil_LowArousal_OffTStamp ] =changepoints(-pupil_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
% title('PupilLowArousal');saveas(h2,fullfile(indivFigureFolder,'PupilLowArousal'));
%get z thresholds based on face data during quiescence only
b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
face_quiescence=face_Norm(pupilTime_quiescence);
zthres_High=quantile(face_quiescence,0.60);
zthres_Low=quantile(face_quiescence,0.40);

%get on and off timestamps for high and low face movment
[h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
title('FaceHighArousal');saveas(h3,fullfile(indivFigureFolder,'FaceHighArousal'));
[h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
title('FaceLowArousal');saveas(h4,fullfile(indivFigureFolder,'FaceLowArousal'));


% get pupil on and face on times if both on and off times occur entirely during sustained queiscence states identified in the previous step
% toDelete=ones(1,length(Pupil_HighArousal_OnTStamp));
% for rj=1:length(Pupil_HighArousal_OnTStamp)
%     tmp = find (Pupil_HighArousal_OnTStamp(rj)>=sitOn_final & Pupil_HighArousal_OffTStamp(rj)<=sitOff_final);
%     toDelete(rj)=isempty(tmp);
% end
% Pupil_HighArousal_On_int=Pupil_HighArousal_OnTStamp(~toDelete)+params.TimeSinceArousalOn;
% Pupil_HighArousal_Off_int=Pupil_HighArousal_OffTStamp(~toDelete)+params.TimeBeforeArousalOff;

% toDelete=ones(1,length(Pupil_LowArousal_OnTStamp));
% for rj=1:length(Pupil_LowArousal_OnTStamp)
%     tmp = find (Pupil_LowArousal_OnTStamp(rj)>=sitOn_final & Pupil_LowArousal_OffTStamp(rj)<=sitOff_final);
%     toDelete(rj)=isempty(tmp);
% end
% Pupil_LowArousal_On_int=Pupil_LowArousal_OnTStamp(~toDelete)+params.TimeSinceArousalOn;
% Pupil_LowArousal_Off_int=Pupil_LowArousal_OffTStamp(~toDelete)+params.TimeBeforeArousalOff;

toDelete=ones(1,length(Face_HighArousal_OnTStamp));
for rj=1:length(Face_HighArousal_OnTStamp)
    tmp = find (Face_HighArousal_OnTStamp(rj)>=sitOn_final & Face_HighArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
Face_HighArousal_On_int=Face_HighArousal_OnTStamp(~toDelete)+params.TimeSinceArousalOn;
Face_HighArousal_Off_int=Face_HighArousal_OffTStamp(~toDelete)+params.TimeBeforeArousalOff;

toDelete=ones(1,length(Face_LowArousal_OnTStamp));
for rj=1:length(Face_LowArousal_OnTStamp)
    tmp = find (Face_LowArousal_OnTStamp(rj)>=sitOn_final & Face_LowArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
Face_LowArousal_On_int=Face_LowArousal_OnTStamp(~toDelete)+params.TimeSinceArousalOn;
Face_LowArousal_Off_int=Face_LowArousal_OffTStamp(~toDelete)+params.TimeBeforeArousalOff;

%get high face/pupil states that are preceded by low face state
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
Face_HighArousal_On_int1=Face_HighArousal_On_int(logical(finalIndex));
Face_HighArousal_Off_int1=Face_HighArousal_Off_int(logical(finalIndex));

% finalIndex=zeros(1,length(Pupil_HighArousal_On_int));
% for ll=1:length(Pupil_HighArousal_On_int)
%     currOnPoint=Pupil_HighArousal_On_int(ll);
%     tmp=currOnPoint-Pupil_LowArousal_Off_int;
%     idxtmp=find(tmp>=0 &tmp<=0.5); %find states where low arousal ends at less than 0.5s before high arousal onset;
%     lowarousallength=Pupil_LowArousal_Off_int(idxtmp)-Pupil_LowArousal_On_int(idxtmp);
%     if length(lowarousallength)>1 || isempty(lowarousallength), continue, end
%     if lowarousallength>=params.TimeBeforeArousalOn
%         finalIndex(ll)=1;
%         Pupil_HighArousal_On_int(ll)=Pupil_LowArousal_Off_int(idxtmp); %replace aorusal on time with low arousal off time
%     end
% end
% Pupil_HighArousal_On_int1=Pupil_HighArousal_On_int(logical(finalIndex));
% Pupil_HighArousal_Off_int1=Pupil_HighArousal_Off_int(logical(finalIndex));


%only select identified high arousal states are at least minimum criterion s long
minimumDuration = params.TimeSinceArousalOn+params.TimeBeforeArousalOff+params.minArousalDuration;
% idx1=find((Pupil_HighArousal_Off_int1-Pupil_HighArousal_On_int1)>=minimumDuration);
% Pupil_On_final=Pupil_HighArousal_On_int1(idx1); Pupil_Off_final=Pupil_HighArousal_Off_int1(idx1);
idx3=find((Face_HighArousal_Off_int1-Face_HighArousal_On_int1)>=minimumDuration);
Face_On_final=Face_HighArousal_On_int1(idx3); Face_Off_final=Face_HighArousal_Off_int1(idx3);
%
%% plot on and off times for sanity check
%plot locomotion state on off times
h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified state

%plot quiescence state on off times
set(0,'CurrentFigure',h);ax1=subplot(3,1,1);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('QuiescenceState');
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

% %plot pupil and face state on off times
% set(0,'CurrentFigure',h);ax2=subplot(3,1,2);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on;  title('PupilHighState');
% for tt=1:length(Pupil_On_final)
%     plot([Pupil_On_final(tt),Pupil_On_final(tt)], ylimits, 'g');
%     plot([Pupil_Off_final(tt),Pupil_Off_final(tt)], ylimits, 'r');
% end

set(0,'CurrentFigure',h);ax3=subplot(3,1,3);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaeState');
for tt=1:length(Face_On_final)
    plot([Face_On_final(tt),Face_On_final(tt)], ylimits, 'g');
    plot([Face_Off_final(tt),Face_Off_final(tt)], ylimits, 'r');
end

linkaxes([ax1, ax3],'x');
saveas(h,fullfile(indivFigureFolder,'BehavioralStatesFinalOnOffTimes'));

% tStamps.pupilOn=Pupil_On_final; tStamps.pupilOff=Pupil_Off_final;
tStamps.faceOn=Face_On_final; tStamps.faceOff=Face_Off_final;

end
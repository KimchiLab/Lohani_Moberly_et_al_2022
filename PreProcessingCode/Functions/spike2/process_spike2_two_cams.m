function [timing, channels_data] = process_spike2_two_cams(cedpath, outputpath,data_smr_time_stamp_filename, fsspike2,channels,minRunDuration,minSitDuration,ITITime,stimDur,stimITI,pupilSR)
%outputs timing in second and continuous data 
%load SON library
CEDS64LoadLib(cedpath)
channels_num=length(fieldnames(channels)); 
% parameters for extracting wheel motion
sCFG.sPARAM.dbWheelDiameterM=0.1524;%in m 
sCFG.sPARAM.db1WheelVRange=[0 5];%cut off range for signal
sCFG.sPARAM.dbWindowLenSec=2; % set the temoporal resolution of the analysis, this seems optimal so far
sCFG.sPARAM.blDoPlot=1; % indicate whether you want to plot or not
sCFG.sPARAM.blDoPlotDuration=1:500000; %plot duration in sample units
%% load event channels 
data = process_spike2_smr2mat('', outputpath, data_smr_time_stamp_filename, channels_num);
[channels_data,wheelOn,wheelOff,h1] = get_channels_data_from_samples(data, channels, sCFG,fsspike2);
saveas(h1,fullfile(outputpath,'WheelOnOffTimesCheck.fig')); 
timing.allwheelon=wheelOn(:); timing.allwheeloff=wheelOff(:); 
threshold1=0.6;threshold2=0.05;
[timing.stimstart,timing.stimend]= VisDetectOnOff(channels_data.diode,threshold1,threshold2,stimDur,stimITI,fsspike2);%visual stim
[timing.mesostart,timing.mesoend]=squaredetect(channels_data.mesoframe,.5);%camera start, end
timing.mesostart=timing.mesostart/fsspike2; timing.mesoend=timing.mesoend/fsspike2; 
[timing.bluestart,timing.blueend]=squaredetect(channels_data.blue,.5);%blue light start, end
timing.bluestart=timing.bluestart/fsspike2;timing.blueend=timing.blueend/fsspike2; 
[timing.uvstart,timing.uvend]=squaredetect(channels_data.uv,.5);%uv light start, end 
timing.uvstart=timing.uvstart/fsspike2; timing.uvend=timing.uvend/fsspike2; 
[airpuffstart,airpuffend]=squaredetect(channels_data.air_puff,.5);%airpuf/electrical stim start,end
airpuffstart=airpuffstart/fsspike2;airpuffend=airpuffend/fsspike2; 
if ~isempty(airpuffstart)
tmp=find(diff(airpuffstart)>5)+1; %make sure airpuffs/estim are spaced by at least 5s 
tmp1=airpuffstart(tmp);tmp2=airpuffend(tmp); 
airpuffstart=[airpuffstart(1); tmp1];
airpuffend=[airpuffend(1); tmp2];%this airpuff end time is not accurate however but we don't need this time 
end 
timing.allairpuffstart=airpuffstart(:); timing.allairpuffend=airpuffend(:); 
threshold1=0.6;threshold2=0.05;
[timing.pupilcamstart,timing.pupilcamend]= PupilDetectOnOff_V2(channels_data.pupil,threshold1,threshold2,1/pupilSR,fsspike2);

%% finetune locomotion and airpuff/electrical stimulation timing 
%only use wheel on with quiescence period of 10s and at least 5s of
%running and not occuring during other events (visual stim and airpuff)
idx=(wheelOn>timing.mesostart(1)+minSitDuration); 
wheelOn_t1=wheelOn(idx); 
wheelOff_t1=wheelOff(idx); 

idx1=find((wheelOff_t1-wheelOn_t1)>=minRunDuration); 
newWheelOnTimes=wheelOn_t1(2:end); 
newWheelOffTimes=wheelOff_t1(1:end-1);
idx2=(find((newWheelOnTimes-newWheelOffTimes)>=minSitDuration))+1;
idx3=intersect(idx1,idx2);                                

wheelOn_int=wheelOn_t1(1,idx3);
wheelOff_int=wheelOff_t1(1,idx3); 

%find wheel on times when airpuff/stim  is not given 
allEvts=[timing.stimstart; airpuffstart]; 
if ~isempty(allEvts)
    allEvts=sort(allEvts,'ascend');
    for i=1:length(allEvts)
        if i==1 && length(allEvts)==1
             index{i}=[(find(wheelOn_int<(allEvts(i)-minRunDuration))), (find(wheelOn_int>(allEvts(i)+ITITime)))];
        elseif i==1   
            index{i}=[(find(wheelOn_int<(allEvts(i)-minRunDuration))), (find(wheelOn_int>(allEvts(i)+ITITime) & wheelOn_int<(allEvts(i+1)-minRunDuration)))];
        elseif i>1 && i<length(allEvts)
            index{i}=find(wheelOn_int>(allEvts(i)+ITITime) & wheelOn_int<(allEvts(i+1)-minRunDuration));
        elseif i==length(allEvts)
            index{i}=find(wheelOn_int>allEvts(i)+ITITime);
        end
    end
    Indices=cell2mat(index);
    wheelOn_final=wheelOn_int(1,Indices);
    wheelOff_final=wheelOff_int(1,Indices);
else
    wheelOn_final=wheelOn_int;
    wheelOff_final=wheelOff_int;
end
timing.wheelOn=wheelOn_final(:); timing.wheelOff=wheelOff_final(:); %locomotion onset and offset times in seconds


%extract airpuff/estim times not contaminated by locomotion 
airpuffStart_final=airpuffstart; 
airpuffEnd_final=airpuffend; 
if ~isempty(airpuffStart_final)
for i=1:length(airpuffstart)
    isRun=find(wheelOn>(airpuffstart(i)-2) & wheelOn<airpuffstart(i)+ITITime);
    ToDelete(i)=~isempty(isRun);
    
end 
airpuffStart_final(ToDelete)=[]; 
airpuffEnd_final(ToDelete)=[]; 
end  
timing.airpuffstart=airpuffStart_final(:); timing.airpuffend=airpuffEnd_final(:); %airpuff/estim on/off times in seconds 




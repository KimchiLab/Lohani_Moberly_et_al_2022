function[cleanData,cleanIdxEEG]= cleanupEEG(data,zThresh1,zThresh2,EventDurAn,preEventDur)
%% clean up EEG trials 
for i=1:length(data)
trials=data{1,i}; 
cleanData{i}=nan(size(trials,1),size(trials,2)); 
%get spectograms and power in three frequency ranges and express data as
%z-scores, then remove trials with artifacts based on z-score thresholds  
win=0.5; step=0.05;
movingwin=[win step]; % set the moving window params
params.Fs=5000; % sampling frequency
params.fpass=[1 100]; % frequencies of
params.tapers=[3 5]; % tapers
params.trialave=0; % average over trials
[S1,t,f]=mtspecgramc(trials,movingwin,params);
lowf=find(f<4);
highf=find(f>60); 
midf=find(f>15 & f<40); 
Lowfreqpower=squeeze(mean(S1(:,lowf ,:),2)); 
highfreqpower= squeeze(mean(S1(:,highf ,:),2)); 
midfreqpower= squeeze(mean(S1(:,midf ,:),2)); 
% baseLowpower=Lowfreqpower(find(t<preEventDur),:); 
% baseHighpower=highfreqpower(find(t<preEventDur),:); 
% baseMidpower=midfreqpower(find(t<preEventDur),:); 
% meanLowPower=mean(baseLowpower(:)); stdLowPower=std(baseLowpower(:)); 
% meanHighPower=mean(baseHighpower(:)); stdHighPower=std(baseHighpower(:)); 
% meanMidPower=mean(baseMidpower(:)); stdMidPower=std(baseMidpower(:)); 
for r=1:size(trials,2)
normLowFreq=(Lowfreqpower(:,r)-mean(Lowfreqpower(:,r)))/std(Lowfreqpower(:,r));    
normHighFreq=(highfreqpower(:,r)-mean(highfreqpower(:,r)))/std(highfreqpower(:,r));
normMidFreq=(midfreqpower(:,r)-mean(midfreqpower(:,r)))/std(midfreqpower(:,r));

% normLowFreq=(Lowfreqpower(:,r)-meanLowPower)/stdLowPower;
% normHighFreq=(highfreqpower(:,r)-meanHighPower)/stdHighPower;
% normMidFreq=(midfreqpower(:,r)-meanMidPower)/stdMidPower;

%figure;subplot(2,1,1); plot(trials(:,r)); subplot(2,1,2);plot(normLowFreq);hold on; plot(normHighFreq); plot(normMidFreq), legend('Low','High','Mid'); hold off; 
cleanIdxEEG{i}(r)=1; 
Low=normLowFreq(find(t>=(preEventDur-EventDurAn(1)) & t<=EventDurAn(2)+preEventDur)); threshLow=Low>zThresh1; 
High=normHighFreq(find(t>=(preEventDur-EventDurAn(1))& t<=EventDurAn(2)+preEventDur));threshHigh=High>zThresh1;  
Mid=normMidFreq(find(t>=(preEventDur-EventDurAn(1)) & t<=EventDurAn(2)+preEventDur));threshMid=Mid>zThresh1; 
combined=threshLow+threshHigh+threshMid;
if any(combined==3)|| any(Low>zThresh2)||any(High>zThresh2)||any(Mid>zThresh2)
    disp('BadEEG_Removing')
    cleanIdxEEG{i}(r) =0; 
 
else 
    disp('GoodEEG')    
    cleanData{i}(:,r)=trials(:,r); 
end 
%pause; close all; 
end 
end 
end 

function[h3,h4]=plotBehaviorAverage(Comb,eventName,preEventWin,params)
%organize behavior data and plot as mean+/- sem
%this function plots behavior mean+/sem around events
h3=figure; set(gcf,'position',[10,10,400,900]);ax1=subplot(5,1,1); ax2=subplot(5,1,2); ax3=subplot(5,1,3); ax4=subplot(5,1,4);ax9=subplot(5,1,5);
h4=figure; set(gcf,'position',[10,10,400,900]);ax5=subplot(5,1,1); ax6=subplot(5,1,2); ax7=subplot(5,1,3); ax8=subplot(5,1,4);ax10=subplot(5,1,5);
%run speed
%combine across animals
runSpeed=Comb.blue.(eventName).RunSpeed;
tstamps=((-preEventWin+1/(params.fsspike2/params.wheelDownsSampleFactor)):1/(params.fsspike2/params.wheelDownsSampleFactor):params.postEventWin)';
RunSpeed_animals=cell2mat(cellfun(@(x) nanmean(x,2),runSpeed,'UniformOutput',false));
RunSpeed_animals=(RunSpeed_animals)*100;%convert to cm
RunSpeed_animals=downsample(RunSpeed_animals,params.wheelDownsSampleFactor); % downsample
plotBehavior(RunSpeed_animals,tstamps,h3,ax1);title(strcat(eventName,'-RunSpeed-Animals')); ylim([-5 20]); ylabel('cm/s');xlabel('Time(s)');xticks(-preEventWin:1:params.postEventWin);
%combine across trials of example mouse
RunSpeed_extrials=Comb.blue.(eventName).RunSpeed{params.Mouse};
RunSpeed_extrials=(RunSpeed_extrials)*100;%convert to cm
RunSpeed_extrials=downsample(RunSpeed_extrials,params.wheelDownsSampleFactor); % downsample
plotBehavior(RunSpeed_extrials,tstamps,h4,ax5);title(strcat(eventName,'-RunSpeed-ExampleMouseTrials')); ylim([-5 20]); ylabel('cm/s');xlabel('Time(s)');xticks(-preEventWin:1:params.postEventWin);

%pupil
if isfield(Comb.blue.(eventName),'pupilNorm')
    %combine across animals
    pupilNorm=Comb.blue.(eventName).pupilNorm;
    tstamps=((-preEventWin+1/(params.fspupilcam)):1/(params.fspupilcam):params.postEventWin)';
    pupil_animals=cell2mat(cellfun(@(x) nanmean(x,2),pupilNorm,'UniformOutput',false));
    plotBehavior(pupil_animals,tstamps,h3,ax2);title(strcat(eventName,'-pupil-Animals')); ylabel('Z-Score(PupilDiameter)');xlabel('Time(s)');xticks(-preEventWin:1:params.postEventWin);
    
    %combine across trials of example mouse
    pupil_extrials=Comb.blue.(eventName).pupilNorm{params.Mouse};
    plotBehavior(pupil_extrials,tstamps,h4,ax6);title(strcat(eventName,'-pupil-ExampleMouseTrials')); ylabel('Z-Score(PupilDiameter)');xlabel('Time(s)');  xticks(-preEventWin:1:params.postEventWin);
end

%face PC1
if isfield(Comb.blue.(eventName),'facePC1CorrNorm')
    %combine across animals
    faceNorm=Comb.blue.(eventName).facePC1CorrNorm;
    tstamps=((-preEventWin+1/(params.fspupilcam)):1/(params.fspupilcam):params.postEventWin)';
    face_animals=cell2mat(cellfun(@(x) nanmean(x,2),faceNorm,'UniformOutput',false));
    plotBehavior(face_animals,tstamps,h3,ax3);title(strcat(eventName,'-face-Animals')); ylabel('Z-Score(facePixels)');xlabel('Time(s)');xticks(-preEventWin:1:params.postEventWin);
    
    %combine across trials of example mouse
    face_extrials=Comb.blue.(eventName).facePC1CorrNorm{params.Mouse};
    plotBehavior(face_extrials,tstamps,h4,ax7);title(strcat(eventName,'-face-ExampleMouseTrials')); ylabel('Z-Score(facePixels)');xlabel('Time(s)');  xticks(-preEventWin:1:params.postEventWin);
end

%EEG
if isfield(Comb.blue.(eventName),'EEG')
    EEG=Comb.blue.(eventName).EEG;
    % clean up trials, make spectrogram, and calculate bandpower in high and low freqs
    zThresh1=2; zThresh2=5; %Z-score thresholds for removing trials with artifacts
    EventDurAn=[preEventWin params.postEventWin];%determine whether EEG is good or bad only during this time interval around event;
    [cleanEEG,cleanIdxEEG]= cleanupEEG(EEG,zThresh1,zThresh2,EventDurAn,preEventWin);
    movingwin=[params.win params.step]; % set the moving window params
    for z=1:size(cleanEEG,2)
        tmp=EEG{z};
        data=tmp(:,logical(cleanIdxEEG{z}));
        [S1,t,f]=mtspecgramc(data,movingwin,params);
        baseS1=mean(S1,1); stdS1=std(S1,0,1);  %normalize to mean activity
        S1norm{z}=(S1-baseS1)./stdS1;
        tmp1=mean(S1norm{z},3);
        meanS{z}=tmp1;
        lowPower{z}=mean(tmp1(:,f<10),2);
        highPower{z}=mean(tmp1(:,f>30),2);
        
        if z==params.Mouse  %extract data for the example mouse
            lowpower_extrials=squeeze(mean(S1norm{z}(:,f<10,:),2));
            highpower_extrials=squeeze(mean(S1norm{z}(:,f>30,:),2));
        end
        
    end
    %take a mean across all animals
    %extract low and high frequency bandpower
    lowpower_animal=cat(2,lowPower{:});
    timestamps=(t-preEventWin)';
    ab1=plotBehavior(lowpower_animal,timestamps,h3,ax4,'b'); ylim([-1 1]); hold on;
    highpower_animal=cat(2,highPower{:});
    ab2=plotBehavior(highpower_animal,timestamps,h3,ax4,'k'); ylabel('Z-score (Power)'); xticks(-preEventWin:1:params.postEventWin); xlabel('Time (s)');
    legend([ab1,ab2],[{'Low'},{'High'}]); title(strcat(eventName,'-EEG-Animals'));
    
    %plot mean spectrogram
    ave_spec=cellfun(@(x) nanmean(x,3),S1norm,'UniformOutput',false);
    combined_spec=cat(3,ave_spec{:});
    meanCombSpec=nanmean(combined_spec,3);
    set(0, 'CurrentFigure', h3);
    set(gcf,'CurrentAxes',ax9); imagesc(timestamps,f,meanCombSpec'); axis xy; ylabel('Frequency'); xlabel('Time(s)'); caxis([-0.5 0.5]);
    
    %take a mean across all trials of example mouse
    ab1=plotBehavior(lowpower_extrials,timestamps,h4,ax8,'b'); ylim([-1 1]); hold on;
    ab2=plotBehavior(highpower_extrials,timestamps,h4,ax8,'k'); ylabel('Z-score (Power)');xticks(-preEventWin:1:params.postEventWin); xlabel('Time (s)');
    legend([ab1,ab2],[{'Low'},{'High'}]); title(strcat(eventName,'-EEG-ExampleMouseTrials'));
    
    %plot mean spectrogram
    meanCombSpec=nanmean(S1norm{params.Mouse},3);
    set(0, 'CurrentFigure', h4);
    set(gcf,'CurrentAxes',ax10); imagesc(timestamps,f,meanCombSpec'); axis xy; ylabel('Frequency'); xlabel('Time(s)');  caxis([-0.5 0.5]);
    
end
end


function [ab]= plotBehavior(Behavior_animals,tstamps,h,ax,color)
%plot behavior data
if nargin==5, c=color; else, c = 'k'; end
meanPop_Behavior_animals=nanmean(Behavior_animals,2);
numAnimals=size(Behavior_animals,2);
stdPop_Behavior_animals=nanstd(Behavior_animals,0,2);
semPop_Behavior_animals=stdPop_Behavior_animals./sqrt(numAnimals);
upperSEM=meanPop_Behavior_animals+semPop_Behavior_animals;
lowerSEM=meanPop_Behavior_animals-semPop_Behavior_animals;
set(0, 'CurrentFigure', h);
set(gcf,'CurrentAxes',ax)
plot(tstamps,meanPop_Behavior_animals,'k');box off;hold on;
xfill=[tstamps; flipud(tstamps)];
yfill=[lowerSEM;flipud(upperSEM)];
ab=fill(xfill,yfill,c,'EdgeColor','None','facealpha',.5) ;

end

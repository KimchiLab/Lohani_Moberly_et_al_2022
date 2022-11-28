function smrxExtractionPreProcessing(tiffsPath,dataSmrxFile,dataVisStimFile,outputPath,channels,dFoF,dFoF_parcells,R,C,params)
%%written by Hadas Benisty 2019 and Sweyta Lohani 2020
%Does the third step of Meso 1P data processing by extracting smrx timestamps for events
%% parameters
if nargin == 10
    fsspike2 = params.fsspike2;
    fsimaging = params.fsimaging;
    pupilSR=params.pupilSR;
    minRunDuration=params.minRunDuration;
    minSitDuration=params.minSitDuration;
    ITITime=params.ITITime;
    preEventWin=params.preEventWin;
    postEventWin=params.postEventWin;
    visStimAn=params.visStimAn;
    visStimDur=params.visDur;
    visStimITI=params.visITI;
else
    fsspike2 = 5e3;
    fsimaging=10;
    pupilSR=10;
    minRunDuration=5;
    minSitDuration=10;
    ITITime=5;
    preEventWin=2;
    postEventWin=5;
    visStimAn=0; 
    visStimDur=2; 
    visStimITI=5; 
    % for blue-uv settings
    params.signalsExtraction.sigs = 'blueuv';
    % for Rcamp AC setting
    params.signalsExtraction.sigs = 'RCaMP_AC';   
end

%% extract timestamps from smrx file and align imaging data to spike2 data
disp(['Extractign smrx timestamps for' tiffsPath]);
if exist(fullfile(outputPath, 'smrx_signals.mat'), 'file')
    if visStimAn
        load(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx','params');
    else
        load(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps');
    end
else
    display(strcat('loading spike 2 file: ',dataSmrxFile));
    cedpath = '../Functions/spike2';
    [timing, channels_data] = process_spike2_two_cams(cedpath, outputPath,dataSmrxFile, fsspike2,channels,minRunDuration,minSitDuration,ITITime,visStimDur,visStimITI,pupilSR);
    
    % get timestamps for all relevant events
    timing.mesostart=timing.mesostart(1:length(timing.mesoend)); % remove the last mesostart if there is an extra one without a corresponding end
    timestamps.timaging = (timing.mesostart+timing.mesoend)/2;
    
    switch params.signalsExtraction.sigs
        case 'blueuv'
            timestamps.timaging=timestamps.timaging(1:2:end);
        case 'RCaMP_AC'
            timestamps.timaging=timestamps.timaging(1:3:end);
    end
    
    timestamps.timaging=timestamps.timaging(params.deterend.filtLen/2:end);%remove the first filtLen/2 timestamps because there is a filtering artifact and those samples are removed in the detrending step
    timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));%remove excess timestamps if you don't have corresponding frames
    timestamps.visStim=timing.stimstart;
    timestamps.airpuff=timing.airpuffstart(timing.airpuffstart>(timestamps.timaging(1)+ preEventWin) & timing.airpuffstart<(timestamps.timaging(end)-postEventWin));%remove events if they occur outside imaging period
    timestamps.wheelOn=timing.wheelOn(timing.wheelOn>(timestamps.timaging(1)+ minSitDuration) & timing.wheelOn<(timestamps.timaging(end)-minRunDuration));%remove events if they occur outside imaging period
    
    %% if you want to analyze visual stim
    if visStimAn
        visData=xlsread(dataVisStimFile);
        contrasts=visData(:,1);
        
        if numel(contrasts)==numel(timestamps.visStim)%get timestamps for each contrast
            visStimidx= find(timing.stimstart>(timestamps.timaging(1)+preEventWin) & timing.stimstart<(timestamps.timaging(end)-postEventWin));
            timestamps.visStim =timing.stimstart(visStimidx);%remove events if they occur outside of imaging period
            contrasts=contrasts(visStimidx); %remove events if they occur outside of imaging period
            uniqContrasts=unique(contrasts);
            for i=1:length(uniqContrasts)
                contrast_Idx{uniqContrasts(i)}=find(contrasts==uniqContrasts(i));
                timestamps.contrasts{uniqContrasts(i)}=timestamps.visStim(contrast_Idx{uniqContrasts(i)});
            end
        else
            error('error:number of stims in excel and spike2 do not mtatch')
        end
        
        save(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx','params');
    else
        save(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','params');
    end
end

%% if visual stimuli were presented, extract ROI based on pixels that are maximally active during 100% contrast presentation
if visStimAn && ~exist(fullfile(outputPath, 'maxVisualTrace.mat'), 'file')
if ~exist('parcells','var'), load('parcells_updated121519.mat'); parcells=parcells_new; end 
V1parcell=parcells.indicators (:,:,2); %mask for left visual cortex 
V1parcellIdx=find(V1parcell==1); 
visActivity_time=(dFoF.green(V1parcellIdx, :));
%take a mean across trials and find peak activity in 0 to 1s from vis stim in the averaged trace 
eventTS=timestamps.contrasts{1,100};preEventWin=2;postEventWin=1; eventTS=eventTS(:)';
[trial_dff]= TrialTimeArrangeDff(visActivity_time,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
meanBaseline=(nanmean(trial_dff(1:params.fsimaging*preEventWin,:,:),1)); 
stdBaseline=(nanstd(trial_dff(1:params.fsimaging*preEventWin,:,:),0,1)); 
ZScore=(trial_dff-meanBaseline)./stdBaseline;
meanActivity=squeeze(mean(ZScore,2));
visActivity=mean(meanActivity(preEventWin*params.fsimaging:end,:),1);

[Ms,Sortidx] = sort(visActivity,'descend');% Sort asending along time dimension
sortedParcellIdx=V1parcellIdx(Sortidx); 
maxVisSortedIdx=sortedParcellIdx(1:ceil(length(Ms)*0.2));%top 20% of values 
maxVisSortedIdx=sort(maxVisSortedIdx);

%make brain color map with overlay of thresholded pixels and parcell boundary
[brain_trial_dff]= TrialTimeArrangeDff(dFoF.green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
tmp=nanmean(brain_trial_dff(preEventWin*params.fsimaging:end,:,:),2);
meanBrainMap=squeeze(nanmean(tmp,1));
reshapedMap=reshape(meanBrainMap,R,C);
h1=figure;subplot(2,1,1);imagesc(reshapedMap);caxis([-0.005 0.007]); colorbar;
meanBrainMap(maxVisSortedIdx)=-1000; top20ImageV=reshape(meanBrainMap,R,C);
subplot(2,1,2);cmap = [ 1 1 1 ; parula(200) ];imagesc(top20ImageV);colormap(cmap);caxis([-0.005 0.007]);  hold on; 
%overlay parcell boundaries
        E=parcells.indicators(:,:,2);
        [B,~] = bwboundaries(E);hold on;
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1),'b', 'LineWidth', 0.1)
        end
hold off;
saveas(h1,fullfile(outputPath,'top20Pixels_visStim'));
save(fullfile(outputPath, 'maxVisualTrace.mat'),'maxVisSortedIdx');
end 
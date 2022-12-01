%% Figure S3 This script calculate peak of negative flourescence dips around a particular event (such as visual stim for GRABs mice or airpuff for GFP control mice)
%written by Sweyta Lohani 2020
close all; clear all
%% user defined folder inputs
inputFolder='W:\GRABS_Data\Analyzed_SVDMethodPatch14\Raw and UVPix outputs for dual mice and egfp\DualMice';%where pixelwise regressed and uncorrected data are
inputFolder1='W:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice';%where spatial SVD data and smrx data are
MiceAnalyze=[{'grab05\'},{'grab06\'},{'grab07\'},{'grab08\'},{'grab09\'},{'grab10\'}];%grabs/rcamp pixelwise regressed and uncorrected animals 
MiceAnalyze1=[{'grabAM05\imaging with 575 excitation\'},{'grabAM06\imaging with 575 excitation\'},{'grabAM07\imaging with 575 excitation\'},{'grabAM08\imaging with 575 excitation\'},{'grabAM09\imaging with 575 excitation\'},{'grabAM10\imaging with 575 excitation\'}]; %grabsrcamp data, spatial SVD correctted
figuresFolder='F:\Figures\VisualStimDipMethodComparison\DualMice';%where figures and output will be saved
%% user-selected input parameters
params.fsimaging=10;%imaging sampling rate
params.preEventWin=2;%pre event duration
params.postEventWin=5;%post event duration
params.baselineWin=2; %use this period to calculate baseline during within trial normalization
%%
addpath(genpath('./Functions'));
if ~exist(figuresFolder),mkdir(figuresFolder); end
V1Idx=2; %parcell idx for V1 in the left hemisphere
%% for each mouse, load spont and airpuff folders
for animal=1:length(MiceAnalyze)
    mainDir =fullfile(inputFolder,MiceAnalyze{animal});
    mainDir1 =fullfile(inputFolder1,MiceAnalyze1{animal});
    folders=dir(mainDir);
    dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
    DirFolders= folders(dirFlags);
    noDrugFlag=~contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true) ;%ignore drug sessions
    DirFolders=DirFolders(noDrugFlag);
    
    folders=dir(mainDir1);
    dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
    DirFolders1= folders(dirFlags);
    noDrugFlag=~contains({DirFolders1.name},'postDrug','IgnoreCase',true)&~contains({DirFolders1.name},'Extra','IgnoreCase',true) ;%ignore drug sessions
    DirFolders1=DirFolders1(noDrugFlag);
    
    %for each subfolder
    for folder =1:length(DirFolders)
        if ~strcmp(DirFolders(folder).name,DirFolders1(folder).name)
            warning('Directory Folders do not match');
            disp(DirFolders(folder).name);
            continue;
        end
        currfolder=fullfile(mainDir, DirFolders(folder).name);
        currfolder1=fullfile(mainDir1, DirFolders1(folder).name);
        % load spike2 and pupil/face data
        load(fullfile(currfolder1,'smrx_signals.mat'),'channels_data','timestamps','timing')
        %load imaging time series (Spatial SVD corrected)
        load(fullfile(currfolder1,'final_dFoF_parcels.mat'),'dFoF_parcells');
        %load imaging time series (pixelwise regression corrected)
        if ~exist(fullfile(currfolder,'dFOF_PixUV_parcells.mat'),'file'), warning('skipping folder'), continue, end;
        load(fullfile(currfolder,'dFOF_PixUV_parcells.mat'),'dFoF_PixUV_parcells');
        %load imaging time series (uncorrected)
        if ~exist(fullfile(currfolder,'dFOF_Raw_parcells.mat'),'file'),warning('skipping folder'), continue, end;
        load(fullfile(currfolder,'dFOF_Raw_parcells.mat'),'dFoF_Raw_parcells');
        
        names=fieldnames(dFoF_PixUV_parcells);
        %for each color
        for rr=1:length(names)
            %extract data from V1 in the blue channel only
            V1data.svdUV=dFoF_parcells.(names{rr})(V1Idx,:);
            V1data.pixUV=dFoF_PixUV_parcells.(names{rr})(V1Idx,:);
            V1data.uncorrected=dFoF_Raw_parcells.(names{rr})(V1Idx,:);
            
            %get trial data for each visual stim (100% contrast)and difference of DFF
            if ~isempty(timestamps.visStim)
                eventTS=timestamps.contrasts{1,100}; eventTS=eventTS(:)';
                
                [V1trial,~]= TrialTimeArrangeDff(V1data.svdUV,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(V1trial(1:params.baselineWin*params.fsimaging,:));
                V1trials.vis.svdUV.(names{rr}){folder}=(V1trial-preEventMean).*100; %mulitply by 100 to get percentages
                
                [V1trial,~]= TrialTimeArrangeDff(V1data.pixUV,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(V1trial(1:params.baselineWin*params.fsimaging,:));
                V1trials.vis.pixUV.(names{rr}){folder}=(V1trial-preEventMean).*100;
                
                [V1trial,~]= TrialTimeArrangeDff(V1data.uncorrected,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(V1trial(1:params.baselineWin*params.fsimaging,:));
                V1trials.vis.uncorrected.(names{rr}){folder}=(V1trial-preEventMean).*100;
            end
            
            %get trial data for each airpuff and Z-score data
            if ~isempty(timestamps.airpuff)
                eventTS=timestamps.airpuff; eventTS=eventTS(:)';
                
                [V1trial,~]= TrialTimeArrangeDff(V1data.svdUV,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(V1trial(1:params.baselineWin*params.fsimaging,:));
                V1trials.airpuff.svdUV.(names{rr}){folder}=(V1trial-preEventMean).*100;
                
                [V1trial,~]= TrialTimeArrangeDff(V1data.pixUV,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(V1trial(1:params.baselineWin*params.fsimaging,:));
                V1trials.airpuff.pixUV.(names{rr}){folder}=(V1trial-preEventMean).*100;
                
                [V1trial,~]= TrialTimeArrangeDff(V1data.uncorrected,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(V1trial(1:params.baselineWin*params.fsimaging,:));
                V1trials.airpuff.uncorrected.(names{rr}){folder}=(V1trial-preEventMean).*100;
            end
        end
    end
    % pool data across sessions for each animal, take an average, and get the peak of negative dip in flourescence
    
    for rr=1:length(names)% for each color
        %% for visual stim
        if isfield(V1trials,'vis')
            combinedData=cell2mat(V1trials.vis.svdUV.(names{rr})); %combine across session
            combinedAverage.vis.svdUV.(names{rr}){animal}=nanmean(combinedData,2);
            postStimAverage=combinedAverage.vis.svdUV.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging);
            minAveFlour.vis.svdUV.(names{rr}){animal}=min(postStimAverage);
            
            combinedData=cell2mat(V1trials.vis.pixUV.(names{rr})); %combine across session
            combinedAverage.vis.pixUV.(names{rr}){animal}=nanmean(combinedData,2);
            postStimAverage=combinedAverage.vis.pixUV.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging);
            minAveFlour.vis.pixUV.(names{rr}){animal}=min(postStimAverage);
            
            combinedData=cell2mat(V1trials.vis.uncorrected.(names{rr})); %combine across session
            combinedAverage.vis.uncorrected.(names{rr}){animal}=nanmean(combinedData,2);
            postStimAverage=combinedAverage.vis.uncorrected.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging);
            minAveFlour.vis.uncorrected.(names{rr}){animal}=min(postStimAverage);
        end
        %% for airpuffs
        if isfield(V1trials,'airpuff')
            combinedData=cell2mat(V1trials.airpuff.svdUV.(names{rr})); %combine across session
            combinedAverage.airpuff.svdUV.(names{rr}){animal}=nanmean(combinedData,2);
            postStimAverage=combinedAverage.airpuff.svdUV.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging);
            minAveFlour.airpuff.svdUV.(names{rr}){animal}=min(postStimAverage);
            
            combinedData=cell2mat(V1trials.airpuff.pixUV.(names{rr})); %combine across session
            combinedAverage.airpuff.pixUV.(names{rr}){animal}=nanmean(combinedData,2);
            postStimAverage=combinedAverage.airpuff.pixUV.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging);
            minAveFlour.airpuff.pixUV.(names{rr}){animal}=min(postStimAverage);
            
            combinedData=cell2mat(V1trials.airpuff.uncorrected.(names{rr})); %combine across session
            combinedAverage.airpuff.uncorrected.(names{rr}){animal}=nanmean(combinedData,2);
            postStimAverage=combinedAverage.airpuff.uncorrected.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging);
            minAveFlour.airpuff.uncorrected.(names{rr}){animal}=min(postStimAverage);
        end
    end
end
% organize data and make paired plots for visual stim
if isfield(minAveFlour,'vis')
    for rr=1:length(names)% for each color
        finalData.vis.data.(names{rr})=[cell2mat(minAveFlour.vis.uncorrected.(names{rr}))',cell2mat(minAveFlour.vis.pixUV.(names{rr}))',cell2mat(minAveFlour.vis.svdUV.(names{rr}))'];
        figure1=figure;
        plot(finalData.vis.data.(names{rr})','-ko','MarkerFaceColor','k'); hold on;
        currMean=mean(finalData.vis.data.(names{rr}),1);
        plot(currMean,'+','MarkerSize',30);  %draw markers at median
        xlim([0.7 3]);xticks(1:3); xticklabels([{'Uncorrected'},{'Pixelwise-UV'},{'SpatialRegres-UV'}]); xtickangle(90); hold on;
        box off; ylabel('MinimumDiffDFF');
        set(gca,'TickDir','out');title(strcat('NegativeFlourescenceDips-VisualStim-',names{rr}));
        saveas(figure1, fullfile(figuresFolder,strcat('NegativeFlourescenceDips-VisualStimPairedPlots-',names{rr})));
        saveas(figure1, fullfile(figuresFolder,strcat('NegativeFlourescenceDips-VisualStimPairedPlots-',names{rr},'.pdf')));
        
        %do one way repeated measures anova followed by posthoc t-tests
        t = table(finalData.vis.data.(names{rr})(:,1),finalData.vis.data.(names{rr})(:,2),finalData.vis.data.(names{rr})(:,3),...
            'VariableNames',{'Method1','Method2','Method3'});
        Method = table([1 2 3]','VariableNames',{'Measurements'});
        RM = fitrm(t,'Method1-Method3~1','WithinDesign',Method);
        finalData.vis.data.stats.(names{rr}).RM= ranova(RM);
        finalData.vis.data.stats.(names{rr}).posthocPVal=nan(1,size(finalData.vis.data.(names{rr}),2));
         finalData.vis.data.stats.(names{rr}).posthocstats=cell(1,size(finalData.vis.data.(names{rr}),2));
        [~,finalData.vis.data.stats.(names{rr}).posthocPVal(1),~,finalData.vis.data.stats.(names{rr}).posthocstats{1}]=ttest(finalData.vis.data.(names{rr})(:,1),finalData.vis.data.(names{rr})(:,2));
        [~,finalData.vis.data.stats.(names{rr}).posthocPVal(2),~,finalData.vis.data.stats.(names{rr}).posthocstats{2}]=ttest(finalData.vis.data.(names{rr})(:,2),finalData.vis.data.(names{rr})(:,3));
        [~,finalData.vis.data.stats.(names{rr}).posthocPVal(3),~,finalData.vis.data.stats.(names{rr}).posthocstats{3}]=ttest(finalData.vis.data.(names{rr})(:,1),finalData.vis.data.(names{rr})(:,3));
    end
end

% organize data and make paired plots for airpuff
if isfield(minAveFlour,'airpuff')
    for rr=1:length(names)% for each color
        finalData.airpuff.data.(names{rr})=[cell2mat(minAveFlour.airpuff.uncorrected.(names{rr}))',cell2mat(minAveFlour.airpuff.pixUV.(names{rr}))',cell2mat(minAveFlour.airpuff.svdUV.(names{rr}))'];
        figure1=figure;
        plot(finalData.airpuff.data.(names{rr})','-ko','MarkerFaceColor','k'); hold on;
        currMean=mean(finalData.airpuff.data.(names{rr}),1);
        plot(currMean,'+','MarkerSize',30);  %draw markers at median
        xlim([0.7 3]);xticks(1:3); xticklabels([{'Uncorrected'},{'Pixelwise-UV'},{'SpatialRegres-UV'}]); xtickangle(90); hold on;
        box off; ylabel('MinimumDiffDFF');
        set(gca,'TickDir','out');title(strcat('NegativeFlourescenceDips-Airpuff-',names{rr}));
        saveas(figure1, fullfile(figuresFolder,strcat('NegativeFlourescenceDips-AirpuffPairedPlots-',names{rr})));
        saveas(figure1, fullfile(figuresFolder,strcat('NegativeFlourescenceDips-AirpuffPairedPlots-',names{rr},'.pdf')));
        
        %do one way repeated measures anova followed by posthoc t-tests
        t = table(finalData.airpuff.data.(names{rr})(:,1),finalData.airpuff.data.(names{rr})(:,2),finalData.airpuff.data.(names{rr})(:,3),...
            'VariableNames',{'Method1','Method2','Method3'});
        Method = table([1 2 3]','VariableNames',{'Measurements'});
        RM = fitrm(t,'Method1-Method3~1','WithinDesign',Method);
        finalData.airpuff.data.stats.(names{rr}).RM= ranova(RM);
        finalData.airpuff.data.stats.(names{rr}).posthocPVal=nan(1,size(finalData.airpuff.data.(names{rr}),2));
        finalData.airpuff.data.stats.(names{rr}).posthocstats=cell(1,size(finalData.airpuff.data.(names{rr}),2));
        [~,finalData.airpuff.data.stats.(names{rr}).posthocPVal(1),~,finalData.airpuff.data.stats.(names{rr}).posthocstats{1}]=ttest(finalData.airpuff.data.(names{rr})(:,1),finalData.airpuff.data.(names{rr})(:,2));
        [~,finalData.airpuff.data.stats.(names{rr}).posthocPVal(2),~,finalData.airpuff.data.stats.(names{rr}).posthocstats{2}]=ttest(finalData.airpuff.data.(names{rr})(:,2),finalData.airpuff.data.(names{rr})(:,3));
        [~,finalData.airpuff.data.stats.(names{rr}).posthocPVal(3),~,finalData.airpuff.data.stats.(names{rr}).posthocstats{3}]=ttest(finalData.airpuff.data.(names{rr})(:,1),finalData.airpuff.data.(names{rr})(:,3));
        
    end
end
%save summary data
save(fullfile(figuresFolder,'SummaryData'),'finalData');

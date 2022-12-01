%% Figure S3, hemodynamics. This script calculates peak of negative flourescence dips around a particular event (such as visual stim for GFP/mcherry contro mice)
%written by Sweyta Lohani 2020
close all; clear all
%% user defined folder inputs
figuresFolder='F:\Figures\VisualStimDipMethodComparison\GFPMcherryMice';%where figures and output will be saved
inputFolder='W:\GRABS_Data\Analyzed_SVDMethodPatch14\GFPMcherryControl';%where pixelwise regressed and uncorrected data are%where pixelwise regressed and uncorrected data are
MiceAnalyze=[{'SLgfp01'},{'SLgfp02'},{'SLgfp03'}]; %mice names
%% user-selected input parameters
params.fsimaging=10;%imaging sampling rate
params.preEventWin=2;
params.postEventWin=5;
params.baselineWin=2; %use this period to calculate baseline during within trial normalization
params.contrast=100; %contrast range from 0 to 100% 
params.fullGrating=0; %grating full (1) or 20 degrees (0)
params.duration=1; % option of 1 or 2s 
params.options='DFF'; %option of doing DFF or ZScore
if strcmp(params.options,'ZScore')
    ylimits=[-3 1];
elseif strcmp(params.options,'DFF')
    ylimits=[-1 0.5]; 
end 
%%
addpath(genpath('./Functions'));
if ~exist(figuresFolder),mkdir(figuresFolder); end
V1Idx=2; S1bIdx=34; M2Idx=52; %parcell idx for V1,S1 and M2 in the left hemisphere
%% for each mouse, load folders 
for animal=1:length(MiceAnalyze)
    mainDir =fullfile(inputFolder,MiceAnalyze{animal});
    folders=dir(mainDir);
    dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
    DirFolders= folders(dirFlags);
    noDrugFlag=~contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true) ;%ignore drug sessions
    DirFolders=DirFolders(noDrugFlag);
    
    %for each subfolder
    for folder =1:length(DirFolders)
        currfolder=fullfile(mainDir, DirFolders(folder).name);
        % load spike2 and pupil/face data
        load(fullfile(currfolder,'smrx_signals.mat'),'channels_data','timestamps','timing','Behavior')
        %load imaging time series (Spatial SVD corrected)
        load(fullfile(currfolder,'final_spatUVdFoF_parcels.mat'));
        if ~exist('dFoF_spatUV_parcells','var')
            dFoF_spatUV_parcells=dFoF_parcells; 
        end 
        %load imaging time series (pixelwise uv regression corrected)
        load(fullfile(currfolder,'final_pixuvdFoF_parcels.mat'),'dFoF_pixuv_parcells');
        %load imaging time series (uncorrected)
        load(fullfile(currfolder,'final_rawdFoF_parcels.mat'),'dFoF_parcells');
        
         %load imaging time series (spatial backscatter regression corrected)
         if exist(fullfile(currfolder,'final_bsdFoF_parcels.mat'),'file') 
         load(fullfile(currfolder,'final_bsdFoF_parcels.mat'),'dFoF_bs_parcells'); 
         end 
         %load imaging time series (pixelwise backscatter regression corrected)
         if exist(fullfile(currfolder,'final_pixbsdFoF_parcels.mat'),'file') 
         load(fullfile(currfolder,'final_pixbsdFoF_parcels.mat'),'dFoF_pixbs_parcells'); 
         end 
        
        names=fieldnames(dFoF_pixuv_parcells);
        %for each color
        for rr=1:length(names)
            %extract data from V1 in the blue channel only
            data.svdUV=dFoF_spatUV_parcells.(names{rr});
            data.pixUV=dFoF_pixuv_parcells.(names{rr});
            data.uncorrected=dFoF_parcells.(names{rr});
            
            %extract visual timestamps for parameters specified above
            stimIdx=find(Behavior.contrasts==params.contrast & Behavior.fullGrating==params.fullGrating & Behavior.duration==params.duration);
            eventTS=timestamps.visStimOn(stimIdx);
            eventTS=eventTS(:)';
            
            %get trial data for each visual stim (100% contrast)and difference of DFF
            if ~isempty(eventTS)  
                [Trials,~]= TrialTimeArrangeDff(data.svdUV,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(Trials(1:params.baselineWin*params.fsimaging,:,:),1);
                preEventStd=nanstd(Trials(1:params.baselineWin*params.fsimaging,:,:),0,1);
                if strcmp(params.options,'DFF')
                trialdata.vis.svdUV.(names{rr}){folder,animal}=(Trials-preEventMean).*100; %mulitply by 100 to get percentages
                elseif strcmp(params.options,'ZScore')
                 trialdata.vis.svdUV.(names{rr}){folder,animal}=(Trials-preEventMean)./preEventStd; 
                end 
                
                [Trials,~]= TrialTimeArrangeDff(data.pixUV,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(Trials(1:params.baselineWin*params.fsimaging,:,:),1);
                preEventStd=nanstd(Trials(1:params.baselineWin*params.fsimaging,:,:),0,1);
                if strcmp(params.options,'DFF')
                    trialdata.vis.pixUV.(names{rr}){folder,animal}=(Trials-preEventMean).*100; %mulitply by 100 to get percentages
                elseif strcmp(params.options,'ZScore')
                    trialdata.vis.pixUV.(names{rr}){folder,animal}=(Trials-preEventMean)./preEventStd; 
                end

                [Trials,~]= TrialTimeArrangeDff(data.uncorrected,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(Trials(1:params.baselineWin*params.fsimaging,:,:),1);
                preEventStd=nanstd(Trials(1:params.baselineWin*params.fsimaging,:,:),0,1);
                if strcmp(params.options,'DFF')
                    trialdata.vis.uncorrected.(names{rr}){folder,animal}=(Trials-preEventMean).*100; %mulitply by 100 to get percentages
                elseif strcmp(params.options,'ZScore')
                   trialdata.vis.uncorrected.(names{rr}){folder,animal}=(Trials-preEventMean)./preEventStd; 
                end
                
                if exist(fullfile(currfolder,'final_bsdFoF_parcels.mat'),'file') 
                data.svdBS=dFoF_bs_parcells.(names{rr});    
                [Trials,~]= TrialTimeArrangeDff(data.svdBS,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(Trials(1:params.baselineWin*params.fsimaging,:,:),1);
                preEventStd=nanstd(Trials(1:params.baselineWin*params.fsimaging,:,:),0,1);
                if strcmp(params.options,'DFF')
                    trialdata.vis.svdBS.(names{rr}){folder,animal}=(Trials-preEventMean).*100; %mulitply by 100 to get percentages
                elseif strcmp(params.options,'ZScore')
                   trialdata.vis.svdBS.(names{rr}){folder,animal}=(Trials-preEventMean)./preEventStd; 
                end
                end 
                
                if exist(fullfile(currfolder,'final_pixbsdFoF_parcels.mat'),'file') 
                data.pixBS=dFoF_pixbs_parcells.(names{rr});
                [Trials,~]= TrialTimeArrangeDff(data.pixBS,timestamps.timaging,params.fsimaging,params.preEventWin,eventTS,params.postEventWin);
                preEventMean=nanmean(Trials(1:params.baselineWin*params.fsimaging,:,:),1);
                preEventStd=nanstd(Trials(1:params.baselineWin*params.fsimaging,:,:),0,1);
                if strcmp(params.options,'DFF')
                    trialdata.vis.pixBS.(names{rr}){folder,animal}=(Trials-preEventMean).*100; %mulitply by 100 to get percentages
                elseif strcmp(params.options,'ZScore')
                   trialdata.vis.pixBS.(names{rr}){folder,animal}=(Trials-preEventMean)./preEventStd; 
                end
                end 
            end
            
          
        end
    end
    % pool data across sessions for each animal, take an average, and get the peak of negative dip in flourescence
    for rr=1:length(names)% for each color
        %% for visual stim
        if isfield(trialdata,'vis')
            combinedData=cell2mat(trialdata.vis.svdUV.(names{rr})(:,animal)'); %combine across session
            combinedAverage.vis.svdUV.(names{rr}){animal}=squeeze(nanmean(combinedData,2));
            postStimData=combinedAverage.vis.svdUV.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging,:);
            minAveFlour.vis.svdUV.(names{rr}){animal}=min(postStimData,[],1);
            minAveFlourV1.vis.svdUV.(names{rr}){animal}=minAveFlour.vis.svdUV.(names{rr}){animal}(V1Idx);
            
            combinedData=cell2mat(trialdata.vis.pixUV.(names{rr})(:,animal)'); %combine across session
            combinedAverage.vis.pixUV.(names{rr}){animal}=squeeze(nanmean(combinedData,2));
            postStimData=combinedAverage.vis.pixUV.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging,:);
            minAveFlour.vis.pixUV.(names{rr}){animal}=min(postStimData,[],1);
            minAveFlourV1.vis.pixUV.(names{rr}){animal}=minAveFlour.vis.pixUV.(names{rr}){animal}(V1Idx);
            
            combinedData=cell2mat(trialdata.vis.uncorrected.(names{rr})(:,animal)'); %combine across session
            combinedAverage.vis.uncorrected.(names{rr}){animal}=squeeze(nanmean(combinedData,2));
            postStimData=combinedAverage.vis.uncorrected.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging,:);
            minAveFlour.vis.uncorrected.(names{rr}){animal}=min(postStimData,[],1);
            minAveFlourV1.vis.uncorrected.(names{rr}){animal}=minAveFlour.vis.uncorrected.(names{rr}){animal}(V1Idx);      
            
            combinedData=cell2mat(trialdata.vis.pixBS.(names{rr})(:,animal)'); %combine across session
            combinedAverage.vis.pixBS.(names{rr}){animal}=squeeze(nanmean(combinedData,2));
            postStimData=combinedAverage.vis.pixBS.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging,:);
            minAveFlour.vis.pixBS.(names{rr}){animal}=min(postStimData,[],1);
            minAveFlourV1.vis.pixBS.(names{rr}){animal}=minAveFlour.vis.pixBS.(names{rr}){animal}(V1Idx);   
            
            combinedData=cell2mat(trialdata.vis.svdBS.(names{rr})(:,animal)'); %combine across session
            combinedAverage.vis.svdBS.(names{rr}){animal}=squeeze(nanmean(combinedData,2));
            postStimData=combinedAverage.vis.svdBS.(names{rr}){animal}(((params.preEventWin*params.fsimaging)+1):(params.postEventWin+params.preEventWin)*params.fsimaging,:);
            minAveFlour.vis.svdBS.(names{rr}){animal}=min(postStimData,[],1);
            minAveFlourV1.vis.svdBS.(names{rr}){animal}=minAveFlour.vis.svdBS.(names{rr}){animal}(V1Idx);   
        end        
    end
end
% organize data and make paired plots for visual stim for V1 and whole brain maps for animal average 
for rr=1:length(names)% for each color
    finalData.vis.data.(names{rr})=[cell2mat(minAveFlourV1.vis.uncorrected.(names{rr}))',cell2mat(minAveFlourV1.vis.pixUV.(names{rr}))',cell2mat(minAveFlourV1.vis.svdUV.(names{rr}))',cell2mat(minAveFlourV1.vis.pixBS.(names{rr}))',cell2mat(minAveFlourV1.vis.svdBS.(names{rr}))'];
    figure1=figure;
    plot(finalData.vis.data.(names{rr})','-ko','MarkerFaceColor','k'); hold on;
    currMean=mean(finalData.vis.data.(names{rr}),1);
    plot(currMean,'+','MarkerSize',30);  %draw markers at median
    xlim([0.7 5]);xticks(1:5); xticklabels([{'Uncorrected'},{'Pixelwise-UV'},{'Spatial-UV'},{'Pixelwise-BS'},{'Spatial-BS'}]); xtickangle(90); hold on;
    box off; ylabel(strcat('Minimum',params.options));
    set(gca,'TickDir','out');title(strcat('NegativeFlourescenceDips-VisualStim-',names{rr}));
    saveas(figure1, fullfile(figuresFolder,strcat('NegativeFlourescenceDips-VisualStimPairedPlots-',names{rr},params.options)));
    saveas(figure1, fullfile(figuresFolder,strcat('NegativeFlourescenceDips-VisualStimPairedPlots-',names{rr},params.options,'.pdf')));
    
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
%% do stats comparing uncorrected negativity between gfp and mcherry and make a paired plot
data=[cell2mat(minAveFlourV1.vis.uncorrected.blue); cell2mat(minAveFlourV1.vis.uncorrected.green)]; 
[~,finalData.vis.data.stats.bluegreenComp.uncorrectedP,~,finalData.vis.data.stats.bluegreenComp.uncorrectedStats]=ttest(data(1,:),data(2,:));
figure2=figure;
plot(data,'-ko','MarkerFaceColor','k'); hold on;
currMean=mean(data,1);
plot(currMean,'+','MarkerSize',30);  %draw markers at median
xlim([0.7 2]);xticks(1:2); xticklabels([{'GFP'},{'Mcherry'}]); xtickangle(90); hold on;
box off; ylabel(strcat('Minimum',params.options));
set(gca,'TickDir','out');title('GFPvsMcherryUncorrectedNegativeDips-VisualStim');
saveas(figure2, fullfile(figuresFolder,strcat('GFPvsMcherryUncorrectedNegativeDips-VisualStim-',params.options)));
saveas(figure2, fullfile(figuresFolder,strcat('GFPvsMcherryUncorrectedNegativeDips-VisualStim-',params.options,'.pdf')));




%% make whole brain negativity plots for gfp and mcherry (uncorrected versus corrected) 
load('parcells_updated121519.mat'); parcells=parcells_new;
Idx.Leftvisual=2:2:16; %visual parcells
Idx.LeftPosteriorParietal=32;% PPC
Idx.LeftRSL_AgL=[18,20];%retrosplenial cortex
Idx.LeftSomat=34:2:48;%somatosensory areas
Idx.LeftFrontalMotor=50:2:52;%frontal-moto area
Idx.Leftauditory=[28,30];%auditory areas

Idx.Rightvisual=1:2:15; %visual parcells
Idx.RightPosteriorParietal=31;% PPC
Idx.RightRSL_AgL=[17,19];%retrosplenial cortex
Idx.RightSomat=33:2:47;%somatosensory areas
Idx.RightFrontalMotor=49:2:51;%frontal-moto area
Idx.Rightauditory=[27,29];%auditory areas
CombinedParcellIdx=[Idx.Leftvisual,Idx.LeftRSL_AgL,Idx.Leftauditory,Idx.LeftPosteriorParietal,Idx.LeftSomat,Idx.LeftFrontalMotor,Idx.Rightvisual,Idx.RightRSL_AgL,Idx.Rightauditory,Idx.RightPosteriorParietal,Idx.RightSomat,Idx.RightFrontalMotor,];
parcellnames=parcells.names(CombinedParcellIdx);
numParcells=length(CombinedParcellIdx); 
V1Idx=1; S1bIdx=14; M2Idx=23; %updated parcell idx for left V1, S1, and M2  

for rr=1:length(names)% for each color
    data1=nanmean(cell2mat(minAveFlour.vis.uncorrected.(names{rr})'),1); 
    data2=nanmean(cell2mat(minAveFlour.vis.pixUV.(names{rr})'),1); 
    data3=nanmean(cell2mat(minAveFlour.vis.svdUV.(names{rr})'),1); 
    data4=nanmean(cell2mat(minAveFlour.vis.pixBS.(names{rr})'),1); 
    data5=nanmean(cell2mat(minAveFlour.vis.svdBS.(names{rr})'),1); 
    [figure1]=plot_parcel_correlation_brainmap(data1(CombinedParcellIdx),'parula',ylimits,CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
    set(gca,'XDir','reverse');title(strcat('Uncorrected-',params.options,names{rr})); 
    [figure2]=plot_parcel_correlation_brainmap(data2(CombinedParcellIdx),'parula',ylimits,CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
    set(gca,'XDir','reverse');title(strcat('PixelwiseUV',params.options,names{rr})); 
    [figure3]=plot_parcel_correlation_brainmap(data3(CombinedParcellIdx),'parula',ylimits,CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
    set(gca,'XDir','reverse');title(strcat('SpatialUV',params.options,names{rr})); 
    [figure4]=plot_parcel_correlation_brainmap(data4(CombinedParcellIdx),'parula',ylimits,CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
    set(gca,'XDir','reverse');title(strcat('PixelwiseBS',params.options,names{rr})); 
    [figure5]=plot_parcel_correlation_brainmap(data5(CombinedParcellIdx),'parula',ylimits,CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
    set(gca,'XDir','reverse');title(strcat('SpatialBS',params.options,names{rr})); 
    
    saveas(figure1, fullfile(figuresFolder,strcat('Uncorrected-',names{rr},params.options)));
    saveas(figure1, fullfile(figuresFolder,strcat('Uncorrected-',names{rr},params.options,'.pdf')));
    saveas(figure2, fullfile(figuresFolder,strcat('PixelwiseUV-',names{rr},params.options)));
    saveas(figure2, fullfile(figuresFolder,strcat('PixelwiseUV-',names{rr},params.options,'.pdf')));
    saveas(figure3, fullfile(figuresFolder,strcat('SpatialUV-',names{rr},params.options)));
    saveas(figure3, fullfile(figuresFolder,strcat('SpatialUV-',names{rr},params.options,'.pdf')));
    saveas(figure4, fullfile(figuresFolder,strcat('PixelwiseBS-',names{rr},params.options)));
    saveas(figure4, fullfile(figuresFolder,strcat('PixelwiseBS-',names{rr},params.options,'.pdf')));
    saveas(figure5, fullfile(figuresFolder,strcat('SpatialBS-',names{rr},params.options)));
    saveas(figure5, fullfile(figuresFolder,strcat('SpatialBS-',names{rr},params.options,'.pdf')));
    
    wholeBrainAveData.(names{rr})=[data1; data2; data3; data4; data5]; 
end 
wholeBrainAveData.Ylabels=parcellnames;
wholeBrainAveData.Xlabels={'uncorrected','pixelwise(395)','Spatial(395)','pixelwise(530)','Spatial(530)'};



%%make trial averaged plots for one example backscatter session inV1
h1=figure; title('V1TrialAverage-ExampleSession')
animal=1; folder=1;
startplot=0;V1Idx=2; 
baselineTime=-(params.preEventWin-1/params.fsimaging):1/params.fsimaging:(params.postEventWin);
for rr=1:length(names)% for each color
    data_svdUV=squeeze(nanmean(trialdata.vis.svdUV.(names{rr}){folder,animal},2));
    data_pixUV=squeeze(nanmean(trialdata.vis.pixUV.(names{rr}){folder,animal},2));
    data_svdBS=squeeze(nanmean(trialdata.vis.svdBS.(names{rr}){folder,animal},2));
    data_pixBS=squeeze(nanmean(trialdata.vis.pixBS.(names{rr}){folder,animal},2));
    data_uncorrected=squeeze(nanmean(trialdata.vis.uncorrected.(names{rr}){folder,animal},2));
   
    
    subplot(2,1,rr);plot(baselineTime,data_uncorrected(:,V1Idx));hold on; plot(baselineTime,data_pixUV(:,V1Idx));plot(baselineTime,data_svdUV(:,V1Idx));
    plot(baselineTime,data_pixBS(:,V1Idx));plot(baselineTime,data_svdBS(:,V1Idx));box off; 
    legend({'Uncorre','PixUV','SpatialUV','PixBS','SpatialBS'});legend box off; 
    title(strcat('V1',names{rr}));  ylabel(params.options);ylim([-1.6 0.7]);
    set(gca,'TickDir','out'); 
end


saveas(h1, fullfile(figuresFolder,'V1TrialAverage-ExampleSession-BackScatterIncluded')); 
saveas(h1, fullfile(figuresFolder,'V1TrialAverage-ExampleSession-BackScatterIncluded.pdf')); 

%save summary data
save(fullfile(figuresFolder,strcat('SummaryData',params.options)),'finalData','wholeBrainAveData','minAveFlour');

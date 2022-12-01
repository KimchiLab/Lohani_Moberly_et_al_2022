%% Figure 3: get amplitudes of sustained staes across parcels and do some statistics for non-pharmacology dual grabs/rcamp mice 
clear all; close all; 
inputFolder='F:\Figures\States';%load sustained state imaging data extracted using the SustainedStatesCorrelatin script 
outputFolder='F:\Figures\States';%where output data will be saved 
if ~exist(outputFolder,'dir'),mkdir(outputFolder), end
load (fullfile(inputFolder,'IndivMouseOutput'),'FaceHighImaging','FaceLowImaging','locImaging')
%parameters
params.signalsExtraction= 'RCaMP_AC'; % 'blueuv (two channels,chat-tigre)' or 'RCaMP_AC (three channels,dual mice)'
params.hem='left'; %choose right or left hemisphere
%% get parcell info 
if strcmp (params.signalsExtraction,'RCaMP_AC')
names={'blue','green'};%channel names 
elseif strcmp (params.signalsExtraction,'blueuv')
names={'blue'}; 
end 
addpath(genpath('F:\Sweyta\GRABS_Data\FinalPaperCodes\ParcellCorrelation'));
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

if strcmp(params.hem,'left')
    parIdx=1:length(CombinedParcellIdx)/2; 
else
    parIdx=length(CombinedParcellIdx)/2+1:length(CombinedParcellIdx);  
end 
leftV1Idx=1;rightV1Idx=24;leftM2Idx=23;rightM2Idx=46;
VisualIdx=[1:8,24:31];%all visual indices left and right hemisphere 
MotorIdx=[22,23,45,46]; %all motor indices left and righ hemisphere 

%% extract mean amplitude for each state 
for rr=1:length(names)
    %% locomotion state
    meanAn.loc.(names{rr})=nan(numParcells,size(locImaging,1)); 
    for jj=1:size(locImaging,1)  
        %get average across epochs within a session 
        for kk=1:size(locImaging,2)
            if ~isempty(locImaging{jj,kk})
                currData=cell2mat(locImaging{jj,kk}.(names{rr}));
                meanSess.loc.(names{rr}){jj,kk}=nanmean(currData,2);               
            end
        end
        %get average across sessions for each animal
        meanAn.loc.(names{rr})(:,jj)=(nanmean(cell2mat(meanSess.loc.(names{rr})(jj,:)),2))*100; %multiply the average df/f amplitude by 100 to get p.c. df/f   
    end
    %% face states
    meanAn.faceH.(names{rr})=nan(numParcells,size(FaceHighImaging,1));   
    meanAn.faceL.(names{rr})=nan(numParcells,size(FaceHighImaging,1));
     for jj=1:size(FaceHighImaging,1)
         %get average across epochs within a session 
        for kk=1:size(FaceHighImaging,2)
            if ~isempty(FaceHighImaging{jj,kk})
                currData=cell2mat(FaceHighImaging{jj,kk}.(names{rr}));
                meanSess.faceH.(names{rr}){jj,kk}=nanmean(currData,2);
                currData=cell2mat(FaceLowImaging{jj,kk}.(names{rr}));
                meanSess.faceL.(names{rr}){jj,kk}=nanmean(currData,2);
            end
        end
        %get average across sessions for each animal
        meanAn.faceH.(names{rr})(:,jj)=(nanmean(cell2mat(meanSess.faceH.(names{rr})(jj,:)),2))*100;   %multiply the average df/f amplitude by 100 to get p.c.df/f       
        meanAn.faceL.(names{rr})(:,jj)=(nanmean(cell2mat(meanSess.faceL.(names{rr})(jj,:)),2))*100;  %multiply the average df/f amplitude by 100 to get  pc. df/f         
    end   
end
%% make mean brain maps for each state- face low, face high and locomotion 
currAve=nanmean(meanAn.loc.blue,2);
[figure1]=plot_parcel_correlation_brainmap(currAve,'parula',[-0.4 0.4],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
title('Locomotion-GRABs');
if strcmp (params.signalsExtraction,'RCaMP_AC')
currAve=nanmean(meanAn.loc.green,2);
[figure4]=plot_parcel_correlation_brainmap(currAve,'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);%brain color map 
title('Locomotion-RCaMP');
end 

currAve=nanmean(meanAn.faceH.blue,2);
[figure2]=plot_parcel_correlation_brainmap(currAve,'parula',[-0.4 0.4],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
title('FaceHigh-GRABs');
if strcmp (params.signalsExtraction,'RCaMP_AC')
currAve=nanmean(meanAn.faceH.green,2);
[figure5]=plot_parcel_correlation_brainmap(currAve,'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);%brain color map 
title('FaceHigh-RCaMP');
end 

currAve=nanmean(meanAn.faceL.blue,2);
[figure3]=plot_parcel_correlation_brainmap(currAve,'parula',[-0.4 0.4],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
title('FaceLow-GRABs');
if strcmp (params.signalsExtraction,'RCaMP_AC')
currAve=nanmean(meanAn.faceL.green,2);
[figure6]=plot_parcel_correlation_brainmap(currAve,'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);%brain color map 
title('FaceLow-RCaMP');
end 

saveas(figure1,fullfile(outputFolder,'BrainMapLoc-GRABs.fig')); saveas(figure1,fullfile(outputFolder,'BrainMapLoc-GRABs.pdf'));
saveas(figure2,fullfile(outputFolder,'BrainMapFaceHigh-GRABs.fig')); saveas(figure2,fullfile(outputFolder,'BrainMapFaceHigh-GRABs.pdf'));
saveas(figure3,fullfile(outputFolder,'BrainMapFaceLow-GRABs.fig')); saveas(figure3,fullfile(outputFolder,'BrainMapFaceLow-GRABs.pdf'));

if strcmp (params.signalsExtraction,'RCaMP_AC')
saveas(figure4,fullfile(outputFolder,'BrainMapLoc-Rcamp.fig')); saveas(figure4,fullfile(outputFolder,'BrainMapLoc-Rcamp.pdf'));
saveas(figure5,fullfile(outputFolder,'BrainMapFaceHigh-Rcamp.fig')); saveas(figure5,fullfile(outputFolder,'BrainMapFaceHigh-Rcamp.pdf'));
saveas(figure6,fullfile(outputFolder,'BrainMapFaceLow-Rcamp.fig')); saveas(figure6,fullfile(outputFolder,'BrainMapFaceLow-Rcamp.pdf'));
end 

%% make brain maps showing average difference between states 
GRABsdata_loc=meanAn.loc.blue-meanAn.faceH.blue; 
currAve=nanmean(GRABsdata_loc,2);
[figure1]=plot_parcel_correlation_brainmap(currAve,'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map 
title('Locomotion-FaceH-GRABs');
if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcampdata_loc=meanAn.loc.green-meanAn.faceH.green; 
currAve=nanmean(Rcampdata_loc,2);
[figure2]=plot_parcel_correlation_brainmap(currAve,'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);%brain color map 
title('Locomotion-FaceH-RCaMP');
end 

GRABsdata_face=meanAn.faceH.blue-meanAn.faceL.blue; 
currAve=nanmean(GRABsdata_face,2);
[figure3]=plot_parcel_correlation_brainmap(currAve,'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);%brain color map 
title('FaceHigh-Low-GRABs');

if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcampdata_face=meanAn.faceH.green-meanAn.faceL.green; 
currAve=nanmean(Rcampdata_face,2);
[figure4]=plot_parcel_correlation_brainmap(currAve,'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);%brain color map 
title('FaceHigh-Low-RCaMP');
end 
saveas(figure1,fullfile(outputFolder,'AverageDifferenceMapLoc-FaceHGRABs.fig')); saveas(figure1,fullfile(outputFolder,'AverageDifferenceMapLoc-FaceHGRABs.pdf'));
saveas(figure3,fullfile(outputFolder,'AverageDifferenceMapFaceHigh-LowGRABs.fig')); saveas(figure3,fullfile(outputFolder,'AverageDifferenceMapFaceHigh-LowGRABs.pdf'));
if strcmp (params.signalsExtraction,'RCaMP_AC')
saveas(figure2,fullfile(outputFolder,'AverageDifferenceMapLoc-FaceHRCaMP.fig')); saveas(figure2,fullfile(outputFolder,'AverageDifferenceMapLoc-FaceHRCaMP.pdf'));
saveas(figure4,fullfile(outputFolder,'AverageDifferenceMapFaceHigh-LowRCaMP.fig')); saveas(figure4,fullfile(outputFolder,'AverageDifferenceMapFaceHigh-LowRCaMP.pdf'));
end 

%% get whole cortex (left hemisphere) averaged value(averaged across parcells) for face high(-face low) and locomotion (-face high) 
GRABs_loc.data=GRABsdata_loc(parIdx,:); GRABs_loc.ave=mean(GRABs_loc.data); %mean across left parcells
GRABs_loc.aveAve=mean(GRABs_loc.ave); GRABs_loc.semAve=std(GRABs_loc.ave)/sqrt(numel(GRABs_loc.ave));%summary stats

GRABs_face.data=GRABsdata_face(parIdx,:); GRABs_face.ave=mean(GRABs_face.data); %mean across left parcells
GRABs_face.aveAve=mean(GRABs_face.ave); GRABs_face.semAve=std(GRABs_face.ave)/sqrt(numel(GRABs_face.ave));%summary stats

if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcamp_loc.data=Rcampdata_loc(parIdx,:); Rcamp_loc.ave=mean(Rcamp_loc.data); %mean across left parcells
Rcamp_loc.aveAve=mean(Rcamp_loc.ave); Rcamp_loc.semAve=std(Rcamp_loc.ave)/sqrt(numel(Rcamp_loc.ave));%summary stats

Rcamp_face.data=Rcampdata_face(parIdx,:); Rcamp_face.ave=mean(Rcamp_face.data); %mean across left parcells
Rcamp_face.aveAve=mean(Rcamp_face.ave); Rcamp_face.semAve=std(Rcamp_face.ave)/sqrt(numel(Rcamp_face.ave));%summary stats
end 

%% get whole cortex (left hemisphere) averaged value(averaged across parcells) for three states without subtraction, face low, face high and locomotion  and do student's t-test for difference from zero
GRABs_locRaw.data=meanAn.loc.blue(parIdx,:); GRABs_locRaw.ave=mean(GRABs_locRaw.data); %mean across left parcells
GRABs_locRaw.aveAve=mean(GRABs_locRaw.ave); GRABs_locRaw.semAve=std(GRABs_locRaw.ave)/sqrt(numel(GRABs_locRaw.ave));%summary stats
[GRABs_locRaw.h, GRABs_locRaw.p, ~, GRABs_locRaw.stats]=ttest(GRABs_locRaw.ave); % one sample t-test

GRABs_faceHRaw.data=meanAn.faceH.blue(parIdx,:); GRABs_faceHRaw.ave=mean(GRABs_faceHRaw.data); %mean across left parcells
GRABs_faceHRaw.aveAve=mean(GRABs_faceHRaw.ave); GRABs_faceHRaw.semAve=std(GRABs_faceHRaw.ave)/sqrt(numel(GRABs_faceHRaw.ave));%summary stats
[GRABs_faceHRaw.h, GRABs_faceHRaw.p,~, GRABs_faceHRaw.stats]=ttest(GRABs_faceHRaw.ave); % one sample t-test

GRABs_faceLRaw.data=meanAn.faceL.blue(parIdx,:); GRABs_faceLRaw.ave=mean(GRABs_faceLRaw.data); %mean across left parcells
GRABs_faceLRaw.aveAve=mean(GRABs_faceLRaw.ave); GRABs_faceLRaw.semAve=std(GRABs_faceLRaw.ave)/sqrt(numel(GRABs_faceLRaw.ave));%summary stats
[GRABs_faceLRaw.h, GRABs_faceLRaw.p,~, GRABs_faceLRaw.stats]=ttest(GRABs_faceLRaw.ave); % one sample t-test

if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcamp_locRaw.data=meanAn.loc.green(parIdx,:); Rcamp_locRaw.ave=mean(Rcamp_locRaw.data); %mean across left parcells
Rcamp_locRaw.aveAve=mean(Rcamp_locRaw.ave); Rcamp_locRaw.semAve=std(Rcamp_locRaw.ave)/sqrt(numel(Rcamp_locRaw.ave));%summary stats
[Rcamp_locRaw.h, Rcamp_locRaw.p,~, Rcamp_locRaw.stats]=ttest(Rcamp_locRaw.ave); % one sample t-test

Rcamp_faceHRaw.data=meanAn.faceH.green(parIdx,:); Rcamp_faceHRaw.ave=mean(Rcamp_faceHRaw.data); %mean across left parcells
Rcamp_faceHRaw.aveAve=mean(Rcamp_faceHRaw.ave); Rcamp_faceHRaw.semAve=std(Rcamp_faceHRaw.ave)/sqrt(numel(Rcamp_faceHRaw.ave));%summary stats
[Rcamp_faceHRaw.h, Rcamp_faceHRaw.p,~, Rcamp_faceHRaw.stats]=ttest(Rcamp_faceHRaw.ave); % one sample t-test

Rcamp_faceLRaw.data=meanAn.faceL.green(parIdx,:); Rcamp_faceLRaw.ave=mean(Rcamp_faceLRaw.data); %mean across left parcells
Rcamp_faceLRaw.aveAve=mean(Rcamp_faceLRaw.ave); Rcamp_faceLRaw.semAve=std(Rcamp_faceLRaw.ave)/sqrt(numel(Rcamp_faceLRaw.ave));%summary stats
[Rcamp_faceLRaw.h, Rcamp_faceLRaw.p,~, Rcamp_faceLRaw.stats]=ttest(Rcamp_faceLRaw.ave); % one sample t-test
end 

%% do two way repeated measures anova with behavioral state (loc, faceh, face low as three levels) and regions as within subjects factors 
numSubjects=size(GRABs_locRaw.data,2); 
numParcells=size(GRABs_locRaw.data,1); 
FACTNAMES={'State','Region'}; 
S=ones(numSubjects*numParcells,1); 
idx=1; 
for i=1:numSubjects
S(idx:numParcells*i,1)=i;
idx=idx+numParcells;
end
S=[S;S;S];

F1=ones(size(S)); 
F1(numParcells*numSubjects+1:(2*numParcells*numSubjects))=2; 
F1((2*numParcells*numSubjects+1):end)=3; 

tmp=1:numParcells; 
F2=repmat(tmp,1,numSubjects*3)'; 

%stats on grabs 
Y=[GRABs_faceLRaw.data(:);GRABs_faceHRaw.data(:);GRABs_locRaw.data(:)]; 
GRABs_RM2.stats = rm_anova2(Y,S,F1,F2,FACTNAMES);

%posthoc on states
[GRABs_RM2.posthoc.faceLfaceH_h, GRABs_RM2.posthoc.faceLfaceH_p, ~,GRABs_RM2.posthoc.faceLfaceH_stats]=ttest(mean(GRABs_faceLRaw.data),mean(GRABs_faceHRaw.data)); 
[GRABs_RM2.posthoc.locfaceH_h, GRABs_RM2.posthoc.locfaceH_p,~,GRABs_RM2.posthoc.locfaceH_stats]=ttest(mean(GRABs_faceHRaw.data),mean(GRABs_locRaw.data)); 

%stats on rcamp 
if strcmp (params.signalsExtraction,'RCaMP_AC')
Y=[Rcamp_faceLRaw.data(:);Rcamp_faceHRaw.data(:);Rcamp_locRaw.data(:)]; 
Rcamp_RM2.stats = rm_anova2(Y,S,F1,F2,FACTNAMES);

%posthoc on states
[Rcamp_RM2.posthoc.faceLfaceH_h, Rcamp_RM2.posthoc.faceLfaceH_p, ~,Rcamp_RM2.posthoc.faceLfaceH_stats]=ttest(mean(Rcamp_faceLRaw.data),mean(Rcamp_faceHRaw.data)); 
[Rcamp_RM2.posthoc.locfaceH_h, Rcamp_RM2.posthoc.locfaceH_p,~,Rcamp_RM2.posthoc.locfaceH_stats]=ttest(mean(Rcamp_faceHRaw.data),mean(Rcamp_locRaw.data)); 
end 

%% compare data between visual and motor areas(Combined left and right hemisphere) for face high,facelow, and locomotion 
startIdx=1; h2=figure; hold on; suptitle('MotorvsVisual'); 
for rr=1:length(names)
%face high
data=meanAn.faceH.(names{rr}); 
an_visual=mean(data(VisualIdx,:));
an_motor=mean(data(MotorIdx,:));

combined_visualmotor.faceH.(names{rr})=[an_motor;an_visual]; 
subplot(2,3,startIdx); title(strcat('facehigh',names{rr}));hold on; 
plot([1,2],combined_visualmotor.faceH.(names{rr}),'.-k','MarkerSize',10);hold on;
currMean=mean(combined_visualmotor.faceH.(names{rr}),2);
plot([1,2],currMean,'+','MarkerSize',30);  %draw markers at mean
xticks(1:2), xticklabels({'Motor','Visual'}); box off; xlim([0.9 2]); ylim([-0.5,1])
set(gca,'Tickdir','out'); 

%do a ttest comparing visual and motor areas 
[~,VisualMotor.faceH.(names{rr}).pval,~,VisualMotor.faceH.(names{rr}).stats]=ttest(an_visual,an_motor); 

%face low
data=meanAn.faceL.(names{rr}); 
an_visual=mean(data(VisualIdx,:));
an_motor=mean(data(MotorIdx,:));

combined_visualmotor.faceL.(names{rr})=[an_motor;an_visual]; 
subplot(2,3,startIdx+1); title(strcat('facelow',names{rr})); hold on; 
plot([1,2],combined_visualmotor.faceL.(names{rr}),'.-k','MarkerSize',10);hold on;
currMean=mean(combined_visualmotor.faceL.(names{rr}),2);
plot([1,2],currMean,'+','MarkerSize',30);  %draw markers at mean
xticks(1:2), xticklabels({'Motor','Visual'}); box off;xlim([0.9 2]); ylim([-0.5,1]) 
set(gca,'Tickdir','out'); 

%do a ttest comparing visual and motor areas 
[~,VisualMotor.faceL.(names{rr}).pval,~,VisualMotor.faceL.(names{rr}).stats]=ttest(an_visual,an_motor); 

%locomotion
data=meanAn.loc.(names{rr}); 
an_visual=mean(data(VisualIdx,:));
an_motor=mean(data(MotorIdx,:));

combined_visualmotor.loc.(names{rr})=[an_motor;an_visual]; 
subplot(2,3,startIdx+2); title(strcat('Locomotion',names{rr}));hold on; 
plot([1,2],combined_visualmotor.loc.(names{rr}),'.-k','MarkerSize',10);hold on;
currMean=mean(combined_visualmotor.loc.(names{rr}),2);
plot([1,2],currMean,'+','MarkerSize',30);  %draw markers at mean
xticks(1:2), xticklabels({'Motor','Visual'}); box off;xlim([0.9 2]); ylim([-0.5,1])
set(gca,'Tickdir','out'); 

%do a ttest comparing visual and motor areas 
[~,VisualMotor.loc.(names{rr}).pval,~,VisualMotor.loc.(names{rr}).stats]=ttest(an_visual,an_motor); 

startIdx=startIdx+3; 
end 

saveas(h2,fullfile(outputFolder,strcat('VisualvsMotorfaceHigh-')));
saveas(h2,fullfile(outputFolder,strcat('VisualvsMotorfaceHigh-','.pdf')));

%% do Spearman's rank correlation on averaged value per parcell with anterior posterior rank 
parcellnames_l=parcellnames(parIdx); [sortedParcellNames,sortIdx]=sort(parcellnames_l); 
%load parcel rank
load('ParcellRank.mat')

%sort data and get a mean across animals 
GRABs_loc.sorteddata=GRABs_loc.data(sortIdx,:); 
GRABs_face.sorteddata=GRABs_face.data(sortIdx,:); 
if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcamp_loc.sorteddata=Rcamp_loc.data(sortIdx,:); 
Rcamp_face.sorteddata=Rcamp_face.data(sortIdx,:); 
end 

grabs_loc=mean(GRABs_loc.sorteddata,2);
grabs_face=mean(GRABs_face.sorteddata,2);
if strcmp (params.signalsExtraction,'RCaMP_AC')
rcamp_loc=mean(Rcamp_loc.sorteddata,2);
rcamp_face=mean(Rcamp_face.sorteddata,2);
end 

%do spearman's correlation 
[GRABs_loc.SpearmanRho,GRABs_loc.SpearmanPVal]=corr(parcellRank,grabs_loc,'Type','Spearman');
[GRABs_face.SpearmanRho,GRABs_face.SpearmanPVal]=corr(parcellRank,grabs_face,'Type','Spearman');

if strcmp (params.signalsExtraction,'RCaMP_AC')
[Rcamp_loc.SpearmanRho,Rcamp_loc.SpearmanPVal]=corr(parcellRank,rcamp_loc,'Type','Spearman');
[Rcamp_face.SpearmanRho,Rcamp_face.SpearmanPVal]=corr(parcellRank,rcamp_face,'Type','Spearman');
end 

%make plots 
[sortedParcellRank,sortIdx]=sort(parcellRank); %sort by anterior posterior rank for plotting 

MeanVal=mean(GRABs_loc.sorteddata,2); SEMVal=(std(GRABs_loc.sorteddata,0,2))./sqrt(size(GRABs_loc.sorteddata,2)); 
SortedMeanVal=MeanVal(sortIdx); SortedSEMVal=SEMVal(sortIdx); 
figure5=figure;subplot(2,2,3);errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('GRABs-loc'); ylim([-0.3 1.3]); box off; 
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out'); 


MeanVal=mean(GRABs_face.sorteddata,2); SEMVal=(std(GRABs_face.sorteddata,0,2))./sqrt(size(GRABs_loc.sorteddata,2)); 
SortedMeanVal=MeanVal(sortIdx); SortedSEMVal=SEMVal(sortIdx); 
subplot(2,2,1);errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('GRABs-face'); ylim([-0.3 1.3]); box off; 
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out'); 

if strcmp (params.signalsExtraction,'RCaMP_AC')
MeanVal=mean(Rcamp_loc.sorteddata,2); SEMVal=(std(Rcamp_loc.sorteddata,0,2))./sqrt(size(Rcamp_loc.sorteddata,2)); 
SortedMeanVal=MeanVal(sortIdx); SortedSEMVal=SEMVal(sortIdx); 
subplot(2,2,4); errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('Rcamp-loc'); ylim([-0.3 1.3]); box off; 
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out'); 


MeanVal=mean(Rcamp_face.sorteddata,2); SEMVal=(std(Rcamp_face.sorteddata,0,2))./sqrt(size(Rcamp_loc.sorteddata,2)); 
SortedMeanVal=MeanVal(sortIdx); SortedSEMVal=SEMVal(sortIdx); 
subplot(2,2,2);errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('Rcamp-face'); ylim([-0.3 1.3]); box off;
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out'); 
end 

saveas(figure5,fullfile(outputFolder,'AnteriorPosteriorGradientSustainedStates.fig')); saveas(figure5,fullfile(outputFolder,'AnteriorPosteriorGradientSustainedStates.pdf'));  

%% compare whole brain average between behavioral states, 
figure6=figure; 
data=[GRABs_face.ave;GRABs_loc.ave];
plot([1,2],data,'-ko','MarkerFaceColor','k'); hold on;
currMean=mean(data,2);
plot([1,2],currMean,'+','MarkerSize',30);  %draw markers at median
xlim([0.75 2.25]);xticks(1:2); xticklabels({'Face','Loc'}); xtickangle(90); hold on;
box off; ylabel('DFF'); set(gca,'TickDir','out');title('Face vs Loc GRABs');
saveas(figure6,fullfile(outputFolder,'PairedPlotsStates-GRABs.fig')); saveas(figure6,fullfile(outputFolder,'PairedPlotsStates-GRABs.pdf'));  

figure8=figure; 
data=[GRABs_faceLRaw.ave;GRABs_faceHRaw.ave; GRABs_locRaw.ave];
plot([1,2,3],data,'-ko','MarkerFaceColor','k'); hold on;
currMean=mean(data,2);
plot([1,2,3],currMean,'+','MarkerSize',30);  %draw markers at median
xlim([0.75 3.25]);xticks(1:3); xticklabels({'FaceL','FaceH','Loc'}); xtickangle(90); hold on;
box off; ylabel('DFF'); set(gca,'TickDir','out');title('FaceL, FaceH vs Loc GRABs');
saveas(figure8,fullfile(outputFolder,'PairedPlotsRawStates-GRABs.fig')); saveas(figure8,fullfile(outputFolder,'PairedPlotsRawStates-GRABs.pdf'));  


if strcmp (params.signalsExtraction,'RCaMP_AC')
figure7=figure; 
data=[Rcamp_face.ave;Rcamp_loc.ave];
plot([1,2],data,'-ko','MarkerFaceColor','k'); hold on;
currMean=median(data,2);
plot([1,2],currMean,'+','MarkerSize',30);  %draw markers at median
xlim([0.75 2.25]);xticks(1:2); xticklabels({'Face','Loc'}); xtickangle(90); hold on;
box off; ylabel('DFF'); set(gca,'TickDir','out');title('Face vs Loc Rcamp');   
saveas(figure7,fullfile(outputFolder,'PairedPlotsStates-Rcamp.fig')); saveas(figure7,fullfile(outputFolder,'PairedPlotsStates-Rcamp.pdf')); 

figure9=figure; 
data=[Rcamp_faceLRaw.ave;Rcamp_faceHRaw.ave; Rcamp_locRaw.ave];
plot([1,2,3],data,'-ko','MarkerFaceColor','k'); hold on;
currMean=mean(data,2);
plot([1,2,3],currMean,'+','MarkerSize',30);  %draw markers at median
xlim([0.75 3.25]);xticks(1:3); xticklabels({'FaceL','FaceH','Loc'}); xtickangle(90); hold on;
box off; ylabel('DFF'); set(gca,'TickDir','out');title('FaceL, FaceH vs Loc GRABs');
saveas(figure9,fullfile(outputFolder,'PairedPlotsRawStates-Rcamp.fig')); saveas(figure9,fullfile(outputFolder,'PairedPlotsRawStates-Rcamp.pdf'));  

end 
%% save output
if strcmp (params.signalsExtraction,'RCaMP_AC')
save(fullfile(outputFolder,'output.mat'),'GRABs_loc','GRABs_face','Rcamp_face','Rcamp_loc','GRABs_locRaw','GRABs_faceLRaw','GRABs_faceHRaw','Rcamp_locRaw','Rcamp_faceLRaw','Rcamp_faceHRaw','GRABs_RM2','Rcamp_RM2','combined_visualmotor','VisualMotor'); 
else
save(fullfile(outputFolder,'output.mat'),'GRABs_loc','GRABs_face','GRABs_locRaw','GRABs_faceLRaw','GRABs_faceHRaw','GRABs_RM2','combined_visualmotor','VisualMotor'); 
end 

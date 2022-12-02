% Figure S7do summary stats and spearman's correlation on transient data (locomotion and face)
clear all; close all; 
inputFolder1='F:\Figures\TransitionStates\Face';%face data folder 
inputFolder2='F:\Figures\TransitionStates\Locomotion';%locomotion data folder 
outputFolder='F:\Figures\TransitionStates'; 
if ~exist(outputFolder,'dir'),mkdir(outputFolder); end 
%parameters
params.signalsExtraction= 'RCaMP_AC'; % 'blueuv' or 'RCaMP_AC'

%% get parcells info 
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
V1Idx=1; S1bIdx=14; M2Idx=23; %updated parcell idx for left V1, S1, and M2  

%% load face data 
load (fullfile(inputFolder1,'summaryData')); 
GRABsdata_face=Blue.face.DiffNorm.allParcells_MaxWin';
GRABsdata_face=GRABsdata_face(CombinedParcellIdx,:); %order with the same indexing as parcellnames 
if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcampdata_face=Green.face.DiffNorm.allParcells_MaxWin';
Rcampdata_face=Rcampdata_face(CombinedParcellIdx,:); %order with the same indexing as parcellnames 
end 

%% load locomotion data 
load (fullfile(inputFolder2,'summaryData')); 
GRABsdata_loc=Blue.wheelOn.DiffNorm.allParcells_MaxWin';
GRABsdata_loc=GRABsdata_loc(CombinedParcellIdx,:); %order with the same indexing as parcellnames 
if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcampdata_loc=Green.wheelOn.DiffNorm.allParcells_MaxWin';
Rcampdata_loc=Rcampdata_loc(CombinedParcellIdx,:); %order with the same indexing as parcellnames 
end 
%% get whole cortex (left hemisphere) averaged value(averaged across parcells) for face and locomotion and do student's t-test for difference from zero
GRABs_loc.data=GRABsdata_loc(1:end/2,:); GRABs_loc.ave=mean(GRABs_loc.data); %mean across left parcells
GRABs_loc.aveAve=mean(GRABs_loc.ave); GRABs_loc.semAve=std(GRABs_loc.ave)/sqrt(numel(GRABs_loc.ave));%summary stats
[GRABs_loc.h, GRABs_loc.p,~,GRABs_loc.stats]=ttest(GRABs_loc.ave); % t-test

GRABs_face.data=GRABsdata_face(1:end/2,:); GRABs_face.ave=mean(GRABs_face.data); %mean across left parcells
GRABs_face.aveAve=mean(GRABs_face.ave); GRABs_face.semAve=std(GRABs_face.ave)/sqrt(numel(GRABs_face.ave));%summary stats
[GRABs_face.h, GRABs_face.p,~,GRABs_face.stats]=ttest(GRABs_face.ave); % t-test

if strcmp (params.signalsExtraction,'RCaMP_AC')
Rcamp_loc.data=Rcampdata_loc(1:end/2,:); Rcamp_loc.ave=mean(Rcamp_loc.data); %mean across left parcells
Rcamp_loc.aveAve=mean(Rcamp_loc.ave); Rcamp_loc.semAve=std(Rcamp_loc.ave)/sqrt(numel(Rcamp_loc.ave));%summary stats
[Rcamp_loc.h, Rcamp_loc.p,~,Rcamp_loc.stats]=ttest(Rcamp_loc.ave); % t-test

Rcamp_face.data=Rcampdata_face(1:end/2,:); Rcamp_face.ave=mean(Rcamp_face.data); %mean across left parcells
Rcamp_face.aveAve=mean(Rcamp_face.ave); Rcamp_face.semAve=std(Rcamp_face.ave)/sqrt(numel(Rcamp_face.ave));%summary stats
[Rcamp_face.h, Rcamp_face.p,~,Rcamp_face.stats]=ttest(Rcamp_face.ave); % t-test
end 

%%  do two way repeated measures ANOVA with behavioral state (locomotion/face as one within subjects factor) and region/parcells as the second within subjects factor 
numSubjects=size(GRABs_loc.data,2); 
numParcells=size(GRABs_loc.data,1); 
FACTNAMES={'State','Region'}; 
S=ones(numSubjects*numParcells,1); 
idx=1; 
for i=1:numSubjects
S(idx:numParcells*i,1)=i;
idx=idx+numParcells;
end
S=[S;S];

F1=ones(size(S)); 
F1(numParcells*numSubjects+1:end)=2; 

tmp=1:numParcells; 
F2=repmat(tmp,1,numSubjects*2)'; 

%stats on grabs 
Y=[GRABs_loc.data(:);GRABs_face.data(:)]; 
GRABs_RM1.stats = rm_anova2(Y,S,F1,F2,FACTNAMES);

%stats on rcamp 
if strcmp (params.signalsExtraction,'RCaMP_AC')
Y=[Rcamp_loc.data(:);Rcamp_face.data(:)]; 
Rcamp_RM1.stats = rm_anova2(Y,S,F1,F2,FACTNAMES);
end 

%% do Spearman's rank correlation on averaged value per parcell with anterior posterior rank 
parcellnames_l=parcellnames(1:end/2); [sortedParcellNames,sortIdx]=sort(parcellnames_l); 
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
figure5=figure;subplot(2,2,3);errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('GRABs-loc'); ylim([0 2.7]); box off; 
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out'); ylabel('DFF'); 


MeanVal=mean(GRABs_face.sorteddata,2); SEMVal=(std(GRABs_face.sorteddata,0,2))./sqrt(size(GRABs_loc.sorteddata,2)); 
SortedMeanVal=MeanVal(sortIdx); SortedSEMVal=SEMVal(sortIdx); 
subplot(2,2,1);errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('GRABs-face'); ylim([0 2.7]); box off; 
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out'); ylabel('DFF'); 

if strcmp (params.signalsExtraction,'RCaMP_AC')
MeanVal=mean(Rcamp_loc.sorteddata,2); SEMVal=(std(Rcamp_loc.sorteddata,0,2))./sqrt(size(Rcamp_loc.sorteddata,2)); 
SortedMeanVal=MeanVal(sortIdx); SortedSEMVal=SEMVal(sortIdx); 
subplot(2,2,4); errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('Rcamp-loc'); ylim([0 2.7]); box off; 
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out');ylabel('DFF');  


MeanVal=mean(Rcamp_face.sorteddata,2); SEMVal=(std(Rcamp_face.sorteddata,0,2))./sqrt(size(Rcamp_loc.sorteddata,2)); 
SortedMeanVal=MeanVal(sortIdx); SortedSEMVal=SEMVal(sortIdx); 
subplot(2,2,2);errorbar(1:23,SortedMeanVal, SortedSEMVal,'Vertical','o','Capsize',0); title('Rcamp-face'); ylim([0 2.7]); box off;
hold on; x=1:23; P = polyfit(x,SortedMeanVal',1);yfit = P(1)*x+P(2); plot(x,yfit,'r');set(gca,'TickDir','out'); ylabel('DFF'); 
end 

saveas(figure5,fullfile(outputFolder,'AnteriorPosteriorGradientTransitionStates.fig')); saveas(figure5,fullfile(outputFolder,'AnteriorPosteriorGradientTransitionStates.pdf'));  

%% save output
if strcmp (params.signalsExtraction,'RCaMP_AC')
save(fullfile(outputFolder,'output.mat'),'GRABs_loc','GRABs_face','Rcamp_face','Rcamp_loc','GRABs_RM1','Rcamp_RM1'); 
else
save(fullfile(outputFolder,'output.mat'),'GRABs_loc','GRABs_face','GRABs_RM1'); 
end 
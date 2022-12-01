%% Figure 3: This script computes amplitudes for sustained states-face high,face low and locomotion for drug and vehicle conditions but focuses on face high only for statistics
close all; clear all
%% user defined inputs
%scopolamine systemic mice
figuresFolder='F:\GRABS_Data\Figures\Scopolamine'; % where figures will be saved
inputFolder='W:\GRABS_Data\ScopolamineAnalyzed';%where input data are located 
MiceAnalyze=[{'CN01'},{'CN02'},{'CN04'},{'HB041_new'},{'HB043'},{'HB044'}];% mice names 
colorChannel={'blue','blue','blue','blue','blue','blue'}; %color channel to extract imaging data in
params.signalsExtraction= 'blueuv'; % 'blueuv' or 'RCaMP_AC'

%% input parameters
params.fsimaging=10;%imaging sampling rate
params.fspupilcam=10; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=10;%for any state, minimum time since any event onset/offset
params.minRunDuration=2;% minimum run duration during locomotion state, previous 5s now 2s
params.minArousalDuration=2; %minimum face/pupil arousal state (high or low arousal),previous 5s now 2s
params.minSitDuration=2;%minimum sit duration during quiescnece state,previous 5s now 2s
params.plotStates=0; %indicate whether a plot of behavioral state timestamps should be generated
params.HemisphereSel='Both'; %indicate whether you want to analyze 'Left' hemisphere only or 'Right' hemisphere only or 'Both' hemispheres
%% add functions
addpath(genpath('./Functions'));
if ~exist(figuresFolder),mkdir(figuresFolder); end

%% get parcell indices for parcells we care about in the left and right hemisphere
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
Idx.RightFrontalMotor=49:2:51;%frontal-motor area
Idx.Rightauditory=[27,29];%auditory areas

if strcmp(params.HemisphereSel,'Left')
    CombinedParcellIdx=[Idx.Leftvisual,Idx.LeftRSL_AgL,Idx.Leftauditory,Idx.LeftPosteriorParietal,Idx.LeftSomat,Idx.LeftFrontalMotor,];
elseif strcmp(params.HemisphereSel,'Right')
    CombinedParcellIdx=[Idx.Rightvisual,Idx.RightRSL_AgL,Idx.Rightauditory,Idx.RightPosteriorParietal,Idx.RightSomat,Idx.RightFrontalMotor,];
elseif strcmp(params.HemisphereSel,'Both')
    CombinedParcellIdx=[Idx.Leftvisual,Idx.LeftRSL_AgL,Idx.Leftauditory,Idx.LeftPosteriorParietal,Idx.LeftSomat,Idx.LeftFrontalMotor,Idx.Rightvisual,Idx.RightRSL_AgL,Idx.Rightauditory,Idx.RightPosteriorParietal,Idx.RightSomat,Idx.RightFrontalMotor,];
end
parcellnames=parcells.names(CombinedParcellIdx);

%specific parcell indices of interest for group analyses 
VisualIdx=[1:8,24:31];%all visual indices left and right hemisphere 
MotorIdx=[22,23,45,46]; %all motor indices left and righ hemisphere 

%indices for brain map 
V1Idx=1;M2Idx=23; S1bIdx=14; 

%% for each mouse, load spont folders corresponding to each drug/vehicle condition and extract state data
if ~exist(fullfile(figuresFolder,strcat('IndivMouseOutput_','Vehicle.mat')),'file')
con='Vehicle'; % whether it's a drug or drug free session
ExtractSustainedState_drugVehicle(con,parcellnames,CombinedParcellIdx,params,MiceAnalyze,inputFolder,figuresFolder);
end 

if ~exist(fullfile(figuresFolder,strcat('IndivMouseOutput_','Drug.mat')),'file')
con='Drug'; %whether it's a drug or drug free session
ExtractSustainedState_drugVehicle(con,parcellnames,CombinedParcellIdx,params,MiceAnalyze,inputFolder,figuresFolder);
end 

%% organize data for each animal 
numParcells=length(CombinedParcellIdx);
vehicle=load(fullfile(figuresFolder,strcat('IndivMouseOutput_','Vehicle')),'locImaging','FaceHighImaging', 'FaceLowImaging');
drug=load(fullfile(figuresFolder,strcat('IndivMouseOutput_','Drug')),'locImaging','FaceHighImaging', 'FaceLowImaging');
veh_loc_final=cell(1,length(MiceAnalyze));veh_faceHigh_final=cell(1,length(MiceAnalyze)); veh_faceLow_final=cell(1,length(MiceAnalyze));
drug_faceLow_final=cell(1,length(MiceAnalyze));drug_faceHigh_final=cell(1,length(MiceAnalyze));drug_loc_final=cell(1,length(MiceAnalyze));
veh_epochs_fH=zeros(1,length(MiceAnalyze)); drug_epochs_fH=zeros(1,length(MiceAnalyze)); 
veh_duration_fH=zeros(1,length(MiceAnalyze));drug_duration_fH=zeros(1,length(MiceAnalyze));
for an=1:length(MiceAnalyze)
    veh_loc_final{an,1}=vehicle.locImaging.Vehicle{an}.(colorChannel{an}); 
    drug_loc_final{an,1}=drug.locImaging.Drug{an}.(colorChannel{an}); 
    veh_faceHigh_final{an,1}=vehicle.FaceHighImaging.Vehicle{an}.(colorChannel{an}); 
    drug_faceHigh_final{an,1}=drug.FaceHighImaging.Drug{an}.(colorChannel{an}); 
    veh_faceLow_final{an,1}=vehicle.FaceLowImaging.Vehicle{an}.(colorChannel{an}); 
    drug_faceLow_final{an,1}=drug.FaceLowImaging.Drug{an}.(colorChannel{an}); 
    
    %get the number of epcohs and duration for face high
    veh_epochs_fH(an)=length(veh_faceHigh_final{an,1});
    drug_epochs_fH(an)=length(drug_faceHigh_final{an,1}); 
    
    veh_duration_fH(an)=sum(cell2mat(cellfun(@(x) size(x,2),veh_faceHigh_final{an,1},'UniformOutput',false))); 
    drug_duration_fH(an)=sum(cell2mat(cellfun(@(x) size(x,2),drug_faceHigh_final{an,1},'UniformOutput',false))); 
end

%% get average dff for each animal for each state for drug and vehicle 
[meanAn_veh.faceH]=averageStateData(veh_faceHigh_final,numParcells);
[meanAn_drug.faceH]=averageStateData(drug_faceHigh_final,numParcells);
[meanAn_veh.faceL]=averageStateData(veh_faceLow_final,numParcells);
[meanAn_drug.faceL]=averageStateData(drug_faceLow_final,numParcells);
[meanAn_veh.loc]=averageStateData(veh_loc_final,numParcells);
[meanAn_drug.loc]=averageStateData(drug_loc_final,numParcells);


%% get mean across animals and make brain maps of average and average difference between drug and vehicle 
h1=figure; suptitle(strcat('Average-Mouse-Parcells'));
subplot(3,5,1);
plot_parcel_correlation_brainmap2(nanmean(meanAn_veh.loc,2),'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('Locomotion-Veh'));

subplot(3,5,2);
plot_parcel_correlation_brainmap2(nanmean(meanAn_veh.faceH,2),'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceHigh-Veh'));

subplot(3,5,3);
plot_parcel_correlation_brainmap2(nanmean(meanAn_veh.faceL,2),'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceLow-Veh'));

% make brain maps showing average difference between states
GRABsdata_loc_pre=meanAn_veh.loc-meanAn_veh.faceH;
subplot(3,5,4);
plot_parcel_correlation_brainmap2(nanmean(GRABsdata_loc_pre,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('Loc-FaceHigh-Veh'))


GRABsdata_face_pre=meanAn_veh.faceH-meanAn_veh.faceL;
subplot(3,5,5);
plot_parcel_correlation_brainmap2(nanmean(GRABsdata_face_pre,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceHigh-Low-Veh'));
%%
subplot(3,5,6);
plot_parcel_correlation_brainmap2(nanmean(meanAn_drug.loc,2),'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('Locomotion-Drug'));

subplot(3,5,7);
plot_parcel_correlation_brainmap2(nanmean(meanAn_drug.faceH,2),'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceHigh-Drug'));

subplot(3,5,8);
plot_parcel_correlation_brainmap2(nanmean(meanAn_drug.faceL,2),'parula',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceLow-Drug'));

% make brain maps showing average difference between states
GRABsdata_loc_post=meanAn_drug.loc-meanAn_drug.faceH;
subplot(3,5,9);
plot_parcel_correlation_brainmap2(nanmean(GRABsdata_loc_post,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('Loc-FaceHigh-Drug'));

GRABsdata_face_post=meanAn_drug.faceH-meanAn_drug.faceL;
subplot(3,5,10);
plot_parcel_correlation_brainmap2(nanmean(GRABsdata_face_post,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceHigh-Low-Drug'));

%%% difference between drug and vehicle
subplot(3,5,11);
plot_parcel_correlation_brainmap2(nanmean(meanAn_drug.loc-meanAn_veh.loc,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('Locomotion-Drug-Vehicle'));

subplot(3,5,12);
plot_parcel_correlation_brainmap2(nanmean(meanAn_drug.faceH-meanAn_veh.faceH,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceHigh-Drug-Vehicle'));

subplot(3,5,13);
plot_parcel_correlation_brainmap2(nanmean(meanAn_drug.faceL-meanAn_veh.faceL,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceLow-Drug-Vehicle'));

% make brain maps showing average difference between states
subplot(3,5,14);
plot_parcel_correlation_brainmap2(nanmean(GRABsdata_loc_post-GRABsdata_loc_pre,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('Loc-FaceHigh-Drug-Vehicle'));

subplot(3,5,15);
plot_parcel_correlation_brainmap2(nanmean(GRABsdata_face_post-GRABsdata_face_pre,2),'redblue',[-0.7 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
set(gca,'XDir','reverse');title(('FaceHigh-Low-Drug-Vehicle'));
saveas(h1,fullfile(figuresFolder,strcat('Average-Mouse-Parcells.fig')));


%% make paired plots for drug versus vehicle for face high only and do statistics 
data_veh=meanAn_veh.faceH; 
data_drug=meanAn_drug.faceH; 
an_vehdrug.all_visual=[nanmean(data_veh(VisualIdx,:)); nanmean(data_drug(VisualIdx,:))]; 
an_vehdrug.all_motor=[nanmean(data_veh(MotorIdx,:)); nanmean(data_drug(MotorIdx,:))]; 

%drug minus vehicle for visual and motor areas 
an_vehsubdrug.all_visual=an_vehdrug.all_visual(2,:)-an_vehdrug.all_visual(1,:); 
an_vehsubdrug.all_motor=an_vehdrug.all_motor(2,:)-an_vehdrug.all_motor(1,:); 

h3=figure; title(strcat('PairedPlots-DrugMinusVehicleVisualvsMotor'));hold on; 
plot([1,2],[an_vehsubdrug.all_motor;an_vehsubdrug.all_visual],'.-k','MarkerSize',10);
currMean=mean([an_vehsubdrug.all_motor;an_vehsubdrug.all_visual],2);
plot([1,2],currMean,'+','MarkerSize',30);  %draw markers at mean
xticks(1:2), xticklabels({'Motor','Visual'}); hold on; box off;
xlim([0.9 2]); ylim([-0.5, 0.2]); 
set(gca,'TickDir','out'); ylabel('DFFDrug-Vehicle'); 
saveas(h3,fullfile(figuresFolder,'DrugminusVehicleVisualvsMotorFaceHigh'));
saveas(h3,fullfile(figuresFolder,'DrugminusVehicleVisualvsMotorFaceHigh.pdf'));

%% do stats with students t-test comparing drug-vehicle against 0 for each region and paired t-test between regions 
[~,stats.visual_P,~,stats.visual_stats]=ttest(an_vehsubdrug.all_visual);
[~,stats.motor_P,~,stats.moto_stats]=ttest(an_vehsubdrug.all_motor);

[~,stats.visualversusmotor_P,~,stats.visualversusmotor_stats]=ttest(an_vehsubdrug.all_motor,an_vehsubdrug.all_visual);


%% do ttest to show that there is no significant difference in number and duration of face high/low states 
[~,stats.duration.pval]=ttest(veh_duration_fH,drug_duration_fH); 
[~,stats.numEpochs.pval]=ttest(veh_epochs_fH,drug_epochs_fH); 
%% save results
save(fullfile(figuresFolder,'SummaryData.mat'),'an_vehdrug','an_vehsubdrug','stats');


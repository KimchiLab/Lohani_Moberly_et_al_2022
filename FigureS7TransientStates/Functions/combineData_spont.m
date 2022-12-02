function [Comb,DataS]= combineData_spont(OverallDir,MiceAnalyze,params,Condition)
%% This function oragnizes imaging data fron spontaneous sessions into trials around locomotion events,
%combines data across sessions for each mouse and pools data across mice
%% inputs
%OverallDir: main input directory containing all the mice subfolders
%MiceAnalyze: mice names
%params: analysis parameters
%Condition: specify if drug was given as 'NoDrug','PreDrug','PostDrug'
%% outputs
%Comb: structure containing outputs for all mice and for all imaging colors, data are combined across sessions
%DataS: structure containing session outputs for all mice and for all imaging colors, contains data from individual sessions
%% imaging color wavelengths
if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
    colors={'green','blue'};
elseif strcmp(params.signalsExtraction.sigs,'blueuv')
    colors={'blue'};
end
%% parcell indices
Idx.Left_Vis_idx=2:2:16; %idx in parcells for left visual areas
Idx.Right_Vis_idx=1:2:15;%right visual areas
Idx.Left_Somat_idx=34:2:48; %idx in parcells for left somatosensory areas
Idx.Right_Somat_idx=33:2:47;%right somatosensory areas
Idx.Left_Frontal_idx=50:2:52; %idx in parcells left for frontal/motor areas (excluded prelimbic and cingulate)
Idx.Right_Frontal_idx=49:2:51;%left frontal/motor areas
Idx.Right_PPC=31;%right posterior parietal cortex
Idx.Left_PPC=32; %left posterior parietal cortex
Idx.RightRSC=[17,19];%right retrosplenial areas
Idx.LeftRSC=[18,20]; %left retrosplenial areas
Idx.LV1=2; %LeftV1 only
Idx.LS1=34; %left somatosensory barrel cortex
Idx.LM2=52;%Left secondary motor area (frontal)
Idx.RV1=1; %Right V1 only
Idx.RS1=33; %Right somatosensory barrel cortex
Idx.RM2=51;%Right secondary motor area (frontal)

%% intitialize final variable that contains combined data for all mice
for colorLen=1:numel(colors)
    currColor=colors{colorLen};
    Comb.(currColor)=[];
    DataS.(currColor)=cell(1,length(MiceAnalyze));
end

%% loop through each mouse's directory
for z=1:length(MiceAnalyze)
    mainDir =fullfile(OverallDir,MiceAnalyze{z});
    folders=dir(mainDir);
    dirFlags = [folders.isdir] & ~contains({folders.name},'vis')& ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');%ignore visual stim sessions
    DirFolders= folders(dirFlags);
    noDrugFlag=~contains({DirFolders.name},'Drug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true);
    vehicleFlag=contains({DirFolders.name},'Vehicle','IgnoreCase',true);
    DrugFlag=contains({DirFolders.name},'Drug','IgnoreCase',true);
    switch Condition
        case 'NoDrug'
            DirFolders=DirFolders(noDrugFlag);
        case 'Vehicle'
            DirFolders=DirFolders(vehicleFlag);
        case 'Drug'
            DirFolders=DirFolders(DrugFlag);
    end
    if isempty(DirFolders),continue, end
    % loop through each session folder and load all data and organize into trials
    for i =1: length(DirFolders)
        currfolder=fullfile(mainDir, DirFolders(i).name);
        %load all imaging and spike2 data
        dFoF=[];
        %load(fullfile(currfolder,'final_dFoF.mat'),'dFoF','R','C'); 
        load(fullfile(currfolder,'final_dFoF_parcels.mat'),'dFoF_parcells');
        load(fullfile(currfolder, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps');
        % load face and pupil data
        files=dir(currfolder);
        fileflags=contains({files.name},'_proc');
        files=files(fileflags);
        if ~isempty(files), load(fullfile(currfolder, files.name),'proc');end
        %extract the current session's data for specified color wavelengths
        for colorLen=1:numel(colors)
            currColor=colors{colorLen};
            %extract locomotion trials
            [DataS.(currColor){z}]=extractSessionData(i,DataS.(currColor){z},currColor,params,dFoF,dFoF_parcells,proc,Idx,timing, channels_data,timestamps.wheelOn,timestamps.timaging,'wheelOn',params.preEventWin,params.postEventWin, params.baselineWin);
        end
        clear dFoF; %clear this variable to reduce memory usage
    end
    %combine data across all sessions for the current mouse
    for colorLen=1:numel(colors)
        currColor=colors{colorLen};
        [Comb.(currColor)]= combineSessionData(z,Comb.(currColor),DataS.(currColor){z},'wheelOn');%wheel on trials
    end
end

%% functions
% first function
    function [DataS] = extractSessionData(sessNum,DataS,color,params,dFoF,dFoF_parcells,proc,Idx,timing, channels_data,evTTimestamps,timaging,eventName,preEventWin,postEventWin, baselineWin)
        % this function extracts data from each session for each color type
        % extract trial delimited activity around events
        disp(strcat('Processing Spontaneous data for',eventName));
        %proces parcellwise imaging data
        eventTS=evTTimestamps(:)';
        [DataS.(eventName).NonNorm_parcellsDFF{sessNum},DataS.(eventName).NonNorm_parcellstimeStamps{sessNum}]= TrialTimeArrangeDff(dFoF_parcells.(color),timaging,params.fsimaging,preEventWin,eventTS,postEventWin);%do parcells only
        DataS.(eventName).NonNorm_parcellsDFF{sessNum}=fillmissing(DataS.(eventName).NonNorm_parcellsDFF{sessNum},'linear',1,'EndValues','nearest');% fill nan values
        
        %process pixelwise imaging data (full video)
%         [DataS.(eventName).NonNorm_pixelsDFF{sessNum},DataS.(eventName).NonNorm_pixelsstimeStamps{sessNum}]= TrialTimeArrangeDff(dFoF.(color),timaging,params.fsimaging,preEventWin,eventTS,postEventWin);%do parcells only
%         DataS.(eventName).NonNorm_pixelsDFF{sessNum}=fillmissing(DataS.(eventName).NonNorm_pixelsDFF{sessNum},'linear',1,'EndValues','nearest');% fill nan values
%         
        %process select parcells
        DataS.(eventName).NonNorm_LeftV1_dffIntp{sessNum}=DataS.(eventName).NonNorm_parcellsDFF{sessNum}(:,:,Idx.LV1);
        DataS.(eventName).NonNorm_RightV1_dffIntp{sessNum}=DataS.(eventName).NonNorm_parcellsDFF{sessNum}(:,:,Idx.RV1);
        DataS.(eventName).NonNorm_LeftS1b_dffIntp{sessNum}=DataS.(eventName).NonNorm_parcellsDFF{sessNum}(:,:,Idx.LS1);
        DataS.(eventName).NonNorm_RightS1b_dffIntp{sessNum}=DataS.(eventName).NonNorm_parcellsDFF{sessNum}(:,:,Idx.RS1);
        DataS.(eventName).NonNorm_LeftM2_dffIntp{sessNum}=DataS.(eventName).NonNorm_parcellsDFF{sessNum}(:,:,Idx.LM2);
        DataS.(eventName).NonNorm_RightM2_dffIntp{sessNum}=DataS.(eventName).NonNorm_parcellsDFF{sessNum}(:,:,Idx.RM2);
        
        baselineFrames=(((baselineWin(1)+preEventWin)*params.fsimaging)+1):((baselineWin(2)+preEventWin)*params.fsimaging); %baseline frames for normalization
        
        %do normalization either as Z-score or as difference from baseline
        baseMean=nanmean(DataS.(eventName).NonNorm_parcellsDFF{sessNum}(baselineFrames,:,:)); baseStd=nanstd(DataS.(eventName).NonNorm_parcellsDFF{sessNum}(baselineFrames,:,:));
        DataS.(eventName).DiffNorm_parcellsDFF{sessNum}= (DataS.(eventName).NonNorm_parcellsDFF{sessNum}-baseMean).*100;%multiply by 100 to get % dff;
        if ~isempty(baseMean)
            DataS.(eventName).ZNorm_parcellsDFF{sessNum}= (DataS.(eventName).NonNorm_parcellsDFF{sessNum}-baseMean)./baseStd;
        else
            DataS.(eventName).ZNorm_parcellsDFF{sessNum}= [];
        end
        
%         baseMean=nanmean(DataS.(eventName).NonNorm_pixelsDFF{sessNum}(baselineFrames,:,:)); baseStd=nanstd(DataS.(eventName).NonNorm_pixelsDFF{sessNum}(baselineFrames,:,:));
%         DataS.(eventName).DiffNorm_pixelsDFF{sessNum}= (DataS.(eventName).NonNorm_pixelsDFF{sessNum}-baseMean).*100;
%         if ~isempty(baseMean)
%             DataS.(eventName).ZNorm_pixelsDFF{sessNum}= (DataS.(eventName).NonNorm_pixelsDFF{sessNum}-baseMean)./baseStd;
%         else
%             DataS.(eventName).ZNorm_pixelsDFF{sessNum}= [];
%         end
%         
        baseMean=nanmean(DataS.(eventName).NonNorm_LeftV1_dffIntp{sessNum}(baselineFrames,:)); baseStd=nanstd(DataS.(eventName).NonNorm_LeftV1_dffIntp{sessNum}(baselineFrames,:));
        DataS.(eventName).DiffNorm_LeftV1_dffIntp{sessNum}= (DataS.(eventName).NonNorm_LeftV1_dffIntp{sessNum}-baseMean).*100;
        if ~isempty(baseMean)
            DataS.(eventName).ZNorm_LeftV1_dffIntp{sessNum}= (DataS.(eventName).NonNorm_LeftV1_dffIntp{sessNum}-baseMean)./baseStd;
        else
            DataS.(eventName).ZNorm_LeftV1_dffIntp{sessNum}= [];
        end
        
        baseMean=nanmean(DataS.(eventName).NonNorm_RightV1_dffIntp{sessNum}(baselineFrames,:)); baseStd=nanstd(DataS.(eventName).NonNorm_RightV1_dffIntp{sessNum}(baselineFrames,:));
        DataS.(eventName).DiffNorm_RightV1_dffIntp{sessNum}= (DataS.(eventName).NonNorm_RightV1_dffIntp{sessNum}-baseMean).*100;
        if ~isempty(baseMean)
            DataS.(eventName).ZNorm_RightV1_dffIntp{sessNum}= (DataS.(eventName).NonNorm_RightV1_dffIntp{sessNum}-baseMean)./baseStd;
        else
            DataS.(eventName).ZNorm_RightV1_dffIntp{sessNum}=[];
        end
        
        baseMean=nanmean(DataS.(eventName).NonNorm_LeftS1b_dffIntp{sessNum}(baselineFrames,:)); baseStd=nanstd(DataS.(eventName).NonNorm_LeftS1b_dffIntp{sessNum}(baselineFrames,:));
        DataS.(eventName).DiffNorm_LeftS1b_dffIntp{sessNum}= (DataS.(eventName).NonNorm_LeftS1b_dffIntp{sessNum}-baseMean).*100;
        if ~isempty(baseMean)
            DataS.(eventName).ZNorm_LeftS1b_dffIntp{sessNum}= (DataS.(eventName).NonNorm_LeftS1b_dffIntp{sessNum}-baseMean)./baseStd;
        else
            DataS.(eventName).ZNorm_LeftS1b_dffIntp{sessNum}=[];
        end
        
        baseMean=nanmean(DataS.(eventName).NonNorm_RightS1b_dffIntp{sessNum}(baselineFrames,:)); baseStd=nanstd(DataS.(eventName).NonNorm_RightS1b_dffIntp{sessNum}(baselineFrames,:));
        DataS.(eventName).DiffNorm_RightS1b_dffIntp{sessNum}= (DataS.(eventName).NonNorm_RightS1b_dffIntp{sessNum}-baseMean).*100;
        if ~isempty(baseMean)
            DataS.(eventName).ZNorm_RightS1b_dffIntp{sessNum}= (DataS.(eventName).NonNorm_RightS1b_dffIntp{sessNum}-baseMean)./baseStd;
        else
            DataS.(eventName).ZNorm_RightS1b_dffIntp{sessNum}=[];
        end
        baseMean=nanmean(DataS.(eventName).NonNorm_LeftM2_dffIntp{sessNum}(baselineFrames,:)); baseStd=nanstd(DataS.(eventName).NonNorm_LeftM2_dffIntp{sessNum}(baselineFrames,:));
        DataS.(eventName).DiffNorm_LeftM2_dffIntp{sessNum}= (DataS.(eventName).NonNorm_LeftM2_dffIntp{sessNum}-baseMean).*100;
        if ~isempty(baseMean)
            DataS.(eventName).ZNorm_LeftM2_dffIntp{sessNum}= (DataS.(eventName).NonNorm_LeftM2_dffIntp{sessNum}-baseMean)./baseStd;
        else
            DataS.(eventName).ZNorm_LeftM2_dffIntp{sessNum}=[];
        end
        
        baseMean=nanmean(DataS.(eventName).NonNorm_RightM2_dffIntp{sessNum}(baselineFrames,:)); baseStd=nanstd(DataS.(eventName).NonNorm_RightM2_dffIntp{sessNum}(baselineFrames,:));
        DataS.(eventName).DiffNorm_RightM2_dffIntp{sessNum}= (DataS.(eventName).NonNorm_RightM2_dffIntp{sessNum}-baseMean).*100;
        if ~isempty(baseMean)
            DataS.(eventName).ZNorm_RightM2_dffIntp{sessNum}= (DataS.(eventName).NonNorm_RightM2_dffIntp{sessNum}-baseMean)./baseStd;
        else
            DataS.(eventName).ZNorm_RightM2_dffIntp{sessNum}=[];
        end
        
        %process behavioral variables
        disp('Processing spike2 data');
        %get run speed
        if isfield(channels_data,'wheelspeed')
            RunSpeed=channels_data.wheelspeed;
            DataS.(eventName).RunSpeed{sessNum}=nan((preEventWin+postEventWin)*params.fsspike2,size(DataS.(eventName).NonNorm_parcellstimeStamps{sessNum},2));
            spike2Time=1/params.fsspike2:1/params.fsspike2:(length(RunSpeed)/params.fsspike2);
            
            for run=1:size(DataS.(eventName).NonNorm_parcellstimeStamps{sessNum},2)
                tmp=RunSpeed(find(spike2Time>=(eventTS(run)-preEventWin) & spike2Time<=(eventTS(run)+postEventWin)));
                DataS.(eventName).RunSpeed{sessNum}(:,run)=tmp(1:size(DataS.(eventName).RunSpeed{sessNum},1));
            end
        end
        
        %extract eeg
        if isfield(channels_data,'EEG')
            EEG=channels_data.EEG;
            DataS.(eventName).EEG{sessNum}=nan((preEventWin+postEventWin)*params.fsspike2,size(DataS.(eventName).NonNorm_parcellstimeStamps{sessNum},2));
            spike2Time=1/params.fsspike2:1/params.fsspike2:(length(EEG)/params.fsspike2);
            
            for ee=1:size(DataS.(eventName).NonNorm_parcellstimeStamps{sessNum},2)
                tmp=EEG(find(spike2Time>=(eventTS(ee)-preEventWin) & spike2Time<=(eventTS(ee)+postEventWin)));
                DataS.(eventName).EEG{sessNum}(:,ee)=tmp(1:size(DataS.(eventName).EEG{sessNum},1));
            end
        end
        
        %extract pupil and face data
        disp('Processing Pupil and face data');
        if exist('proc','var')
            %get normalized pupil/face output if it exists
            if isfield(proc,'output')
                pupilCamTimes=(timing.pupilcamstart+timing.pupilcamend)./2;
                %[DataS.(eventName).pupilNorm{sessNum},~]= TrialTimeArrangeDff(proc.output.pupilNorm(:)',pupilCamTimes,params.fspupilcam,preEventWin,eventTS,postEventWin);
                [DataS.(eventName).facePC1CorrNorm{sessNum},~]= TrialTimeArrangeDff(proc.output.facePC1CorrNorm(:)',pupilCamTimes,params.fspupilcam,preEventWin,eventTS,postEventWin);
                %[DataS.(eventName).whiskerPC1CorrNorm{sessNum},~]= TrialTimeArrangeDff(proc.output.whiskerPC1CorrNorm(:)',pupilCamTimes,params.fspupilcam,preEventWin,eventTS,postEventWin);
            end
        end
        
    end

%% second function
    function[Comb]= combineSessionData(aniNum,Comb,DataS,eventName)
        %this function concatenates data across sessions within each mouse
        Comb.(eventName).ZNorm_parcellsDFF{aniNum}=cell2mat(DataS.(eventName).ZNorm_parcellsDFF);
        Comb.(eventName).ZNorm_LeftV1_dffIntp{aniNum}=cell2mat(DataS.(eventName).ZNorm_LeftV1_dffIntp);
        Comb.(eventName).ZNorm_RightV1_dffIntp{aniNum}=cell2mat(DataS.(eventName).ZNorm_RightV1_dffIntp);
        Comb.(eventName).ZNorm_LeftS1b_dffIntp{aniNum}=cell2mat(DataS.(eventName).ZNorm_LeftS1b_dffIntp);
        Comb.(eventName).ZNorm_RightS1b_dffIntp{aniNum}=cell2mat(DataS.(eventName).ZNorm_RightS1b_dffIntp);
        Comb.(eventName).ZNorm_LeftM2_dffIntp{aniNum}=cell2mat(DataS.(eventName).ZNorm_LeftM2_dffIntp);
        Comb.(eventName).ZNorm_RightM2_dffIntp{aniNum}=cell2mat(DataS.(eventName).ZNorm_RightM2_dffIntp);
        
        Comb.(eventName).NonNorm_parcellsDFF{aniNum}=cell2mat(DataS.(eventName).NonNorm_parcellsDFF);
        Comb.(eventName).NonNorm_LeftV1_dffIntp{aniNum}=cell2mat(DataS.(eventName).NonNorm_LeftV1_dffIntp);
        Comb.(eventName).NonNorm_RightV1_dffIntp{aniNum}=cell2mat(DataS.(eventName).NonNorm_RightV1_dffIntp);
        Comb.(eventName).NonNorm_LeftS1b_dffIntp{aniNum}=cell2mat(DataS.(eventName).NonNorm_LeftS1b_dffIntp);
        Comb.(eventName).NonNorm_RightS1b_dffIntp{aniNum}=cell2mat(DataS.(eventName).NonNorm_RightS1b_dffIntp);
        Comb.(eventName).NonNorm_LeftM2_dffIntp{aniNum}=cell2mat(DataS.(eventName).NonNorm_LeftM2_dffIntp);
        Comb.(eventName).NonNorm_RightM2_dffIntp{aniNum}=cell2mat(DataS.(eventName).NonNorm_RightM2_dffIntp);
        
        Comb.(eventName).DiffNorm_parcellsDFF{aniNum}=cell2mat(DataS.(eventName).DiffNorm_parcellsDFF);
        Comb.(eventName).DiffNorm_LeftV1_dffIntp{aniNum}=cell2mat(DataS.(eventName).DiffNorm_LeftV1_dffIntp);
        Comb.(eventName).DiffNorm_RightV1_dffIntp{aniNum}=cell2mat(DataS.(eventName).DiffNorm_RightV1_dffIntp);
        Comb.(eventName).DiffNorm_LeftS1b_dffIntp{aniNum}=cell2mat(DataS.(eventName).DiffNorm_LeftS1b_dffIntp);
        Comb.(eventName).DiffNorm_RightS1b_dffIntp{aniNum}=cell2mat(DataS.(eventName).DiffNorm_RightS1b_dffIntp);
        Comb.(eventName).DiffNorm_LeftM2_dffIntp{aniNum}=cell2mat(DataS.(eventName).DiffNorm_LeftM2_dffIntp);
        Comb.(eventName).DiffNorm_RightM2_dffIntp{aniNum}=cell2mat(DataS.(eventName).DiffNorm_RightM2_dffIntp);
        
        
        %concatenate data across all sessions and contrasts for behavior
        if isfield(DataS.(eventName),'RunSpeed'),Comb.(eventName).RunSpeed{aniNum}=cell2mat(DataS.(eventName).RunSpeed);end
        if isfield(DataS.(eventName),'EEG'), Comb.(eventName).EEG{aniNum}=cell2mat(DataS.(eventName).EEG);end
        if isfield(DataS.(eventName),'pupilNorm'), Comb.(eventName).pupilNorm{aniNum}=cell2mat(DataS.(eventName).pupilNorm);end
        if isfield(DataS.(eventName),'facePC1CorrNorm'), Comb.(eventName).facePC1CorrNorm{aniNum}=cell2mat(DataS.(eventName).facePC1CorrNorm);end
        if isfield(DataS.(eventName),'whiskerPC1CorrNorm'), Comb.(eventName).whiskerPC1CorrNorm{aniNum}=cell2mat(DataS.(eventName).whiskerPC1CorrNorm);end
        
    end

end


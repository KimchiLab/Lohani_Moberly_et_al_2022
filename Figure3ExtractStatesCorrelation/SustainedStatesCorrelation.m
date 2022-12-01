%% Figures 3,4,5 This script first extracts sustained states (face low, face high, locomotion) and calculate correlations in those states
%(cross-correlation with specified time lag or pearson's correlation) between parcells for different color combinations (blue-blue, green-green)
%or within parcell for two colors (Blue and green). All time are in seconds
%written by Sweyta Lohani 2020
close all; clear all
%% user defined folder inputs and outputs 
figuresFolder='F:\Figures\States'; %where output figures/data will be saved
inputFolder='W:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice';%dual grabs/rcamp mice input data 
MiceAnalyze=[{'grabAM05\imaging with 575 excitation\'},{'grabAM06\imaging with 575 excitation\'},{'grabAM07\imaging with 575 excitation\'},{'grabAM08\imaging with 575 excitation\'},{'grabAM09\imaging with 575 excitation\'},{'grabAM10\imaging with 575 excitation\'}];
Condition='NoDrug'; %'NoDrug','PreDrug','PostDrug' }]; % whethere it's a drug or drug free session
%% user-selected input parameters
params.signalsExtraction= 'RCaMP_AC'; % 'blueuv' or 'RCaMP_AC'
params.fsimaging=10;%imaging sampling rate
params.fspupilcam=10; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=10;%for any state, minimum time since any event onset/offset
params.minRunDuration=5;% minimum run duration during locomotion state
params.minArousalDuration=5; %minimum face/pupil arousal state (high or low arousal)
params.minSitDuration=5;%minimum sit duration during quiescnece state
params.corrType='crosscorr'; %choose between 'pearsons'correlation and time-lagged cross-correlation 'crosscorr'
params.CorrLags=5; %lags for cross-correlation in frames, so approximately +/- 0.5s
params.ShuffleNum=1; %number of times to shuffle data (if you want to use shuffled control for comparison(this is for visual check only, not used in any stats)
params.numPerms=10000; %number of permutations (for permutation stats test)
params.plotStates=1; %indicate whether a plot of behavioral state timestamps should be generated 
params.HemisphereSel='Both'; %indicate whether you want to analyze 'Left' hemisphere only or 'Right' hemisphere only or 'Both' hemispheres  
%% add functions
addpath(genpath('.\Functions'));
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
Idx.RightFrontalMotor=49:2:51;%frontal-moto area
Idx.Rightauditory=[27,29];%auditory areas

if strcmp(params.HemisphereSel,'Left')
CombinedParcellIdx=[Idx.Leftvisual,Idx.LeftRSL_AgL,Idx.Leftauditory,Idx.LeftPosteriorParietal,Idx.LeftSomat,Idx.LeftFrontalMotor,];
elseif strcmp(params.HemisphereSel,'Right')
CombinedParcellIdx=[Idx.Rightvisual,Idx.RightRSL_AgL,Idx.Rightauditory,Idx.RightPosteriorParietal,Idx.RightSomat,Idx.RightFrontalMotor,];   
elseif strcmp(params.HemisphereSel,'Both')
CombinedParcellIdx=[Idx.Leftvisual,Idx.LeftRSL_AgL,Idx.Leftauditory,Idx.LeftPosteriorParietal,Idx.LeftSomat,Idx.LeftFrontalMotor,Idx.Rightvisual,Idx.RightRSL_AgL,Idx.Rightauditory,Idx.RightPosteriorParietal,Idx.RightSomat,Idx.RightFrontalMotor,];
end 
parcellnames=parcells.names(CombinedParcellIdx);
V1Idx=1; S1bIdx=14; M2Idx=23; %updated parcell idx for left V1, S1, and M2  
%% for each mouse, load spont and airpuff folders and perform correlations
for animal=1:length(MiceAnalyze)
    mainDir1 =fullfile(inputFolder,MiceAnalyze{animal});
    folders=dir(mainDir1);
    dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
    DirFolders= folders(dirFlags);
    noDrugFlag=~contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim and drug sessions
    preDrugFlag=contains({DirFolders.name},'preDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
    postDrugFlag=contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
    switch Condition
        case 'NoDrug'
            DirFolders=DirFolders(noDrugFlag);
        case 'PreDrug'
            DirFolders=DirFolders(preDrugFlag);
        case 'PostDrug'
            DirFolders=DirFolders(postDrugFlag);
    end
    
    if isempty(DirFolders),continue, end
    %for each subfolder
    for folder =1:length(DirFolders)
        currfolder=fullfile(mainDir1, DirFolders(folder).name);
        %load imaging time series
        load(fullfile(currfolder,'final_dFoF_parcels.mat'),'dFoF_parcells');
        % load spike2 and pupil/face data
        load(fullfile(currfolder,'smrx_signals.mat'),'channels_data','timestamps','timing')
        files=dir(currfolder);
        for t=1:length(files)
            if contains(files(t).name,'proc','IgnoreCase',true)
                files=files(t);break;
            end
        end
        load(fullfile(currfolder,files.name),'proc')
        face_Norm=proc.output.facePC1CorrNorm;
        wheel_speed = channels_data.wheelspeed;
        
        % get pupil wheel and imaging times
        wheel_time = (1:length(wheel_speed))/params.fsspike2;
        imaging_time = timestamps.timaging;
        pupil_time = timing.pupilcamstart(1:length(face_Norm));
        
        %% extract imaging in selected parcells
        numShuffles=params.ShuffleNum;% number of random surrogate control data to build
        names=fieldnames(dFoF_parcells);
        for i=1:length(names)
            dFoF_parcells.(names{i})=dFoF_parcells.(names{i}) (CombinedParcellIdx,:); % extract data from specific parcells only
            %for each parcell, create shuffled control data by building 1000 surrogate datasets. To do this, data are randomly partionted into two slices and the order of slices are rearranged.
            dFoF_parcells_Shuffled.(names{i})=nan(size(dFoF_parcells.(names{i}),1),size(dFoF_parcells.(names{i}),2),numShuffles);
            for par=1:size(dFoF_parcells.(names{i}),1)
                for numshuff=1:numShuffles
                    data=dFoF_parcells.(names{i})(par,:);
                    frames=100:length(data)-100; %offsets for the splitpoints 
                    splitpoint=frames(randi(length(frames)));
                    tmp1=data(1,1:splitpoint);  tmp2=data(1,splitpoint+1:end);
                    final_data=cat(2,tmp2,tmp1); %break the data into two halves and rearrange the order (circular shifting)
                    dFoF_parcells_Shuffled.(names{i})(par,:,numshuff)=final_data;
                end
            end
        end
        %% make output figures folder
        indivFigureFolder=fullfile(figuresFolder,MiceAnalyze{animal},DirFolders(folder).name);
        if ~exist(indivFigureFolder,'dir'),mkdir(indivFigureFolder), end
        
        %% extract sustained locomotin,face high, and face low states
        [wheelOn_final,wheelOff_final,Face_HighArousal_On_final,Face_HighArousal_Off_final,Face_LowArousal_On_final,Face_LowArousal_Off_final]...
            =getSustainedStates(timing,params,imaging_time,pupil_time,wheel_time, face_Norm,wheel_speed, files,indivFigureFolder);
      %% ensure for comparison purpose, wihitn each session, the number of face low (sit), face high (sit), and locomotion trials are matched. If no wheel or face high or face low trials exist, skip
        if ~isempty(wheelOn_final) &&~isempty(Face_HighArousal_On_final) && ~isempty(Face_LowArousal_On_final)
            %sort epochs by duration and use the same number of epochs across states 
            wheelOnDur=(wheelOff_final-wheelOn_final);
            Face_HighDur=Face_HighArousal_Off_final-Face_HighArousal_On_final;
            Face_LowDur=Face_LowArousal_Off_final-Face_LowArousal_On_final;
            numWheel=length(wheelOn_final); numFaceH=length(Face_HighArousal_On_final); numFaceL=length(Face_LowArousal_On_final);
            totalTrials_wheel=min([numWheel,numFaceH,numFaceL]);
            [sortedWheelDur,sortedWheelidx]=sort(wheelOnDur,'descend');  
            [sortedFaceHighDur,sortedFaceHighidx]=sort(Face_HighDur,'descend');
            [sortedFaceLowDur,sortedFaceLowidx]=sort(Face_LowDur,'descend');
            sortedWheelIdx=sortedWheelidx(1:totalTrials_wheel); 
            sortedFaceHighIdx=sortedFaceHighidx(1:totalTrials_wheel); 
            sortedFaceLowIdx=sortedFaceLowidx(1:totalTrials_wheel); 
            
            wheelOn_final1=wheelOn_final(sortedWheelIdx);
            wheelOff_final1=wheelOff_final(sortedWheelIdx);
            Face_HighArousal_On_final1=Face_HighArousal_On_final(sortedFaceHighIdx);
            Face_HighArousal_Off_final1=Face_HighArousal_Off_final(sortedFaceHighIdx);
            Face_LowArousal_On_final1=Face_LowArousal_On_final(sortedFaceLowIdx);
            Face_LowArousal_Off_final1=Face_LowArousal_Off_final(sortedFaceLowIdx);
            
            %match the duration of epochs across states 
            wheelOnDur=(wheelOff_final1-wheelOn_final1);
            Face_HighArousal_OnDur=Face_HighArousal_Off_final1-Face_HighArousal_On_final1;
            Face_LowArousal_OnDur=Face_LowArousal_Off_final1-Face_LowArousal_On_final1;
            minDur=min([Face_HighArousal_OnDur',Face_LowArousal_OnDur',wheelOnDur'],[],2);
            wheelOff_final1=wheelOn_final1+minDur';
            Face_HighArousal_Off_final1=Face_HighArousal_On_final1+minDur';
            Face_LowArousal_Off_final1=Face_LowArousal_On_final1+minDur';
            
            CompData.loc.TotalDur{animal,folder}=sum(minDur);
            CompData.loc.numEpochs{animal,folder}=totalTrials_wheel;
            
            if strcmp(params.corrType,'crosscorr')
            % run cross correlations between areas on data from three states  and extract the max correlation
            [loccorrMat{animal,folder},loccorrMat_Shuffled{animal,folder},FaceHighcorrMat{animal,folder},FaceHighcorrMat_Shuffled{animal,folder},FaceLowcorrMat{animal,folder},FaceLowcorrMat_Shuffled{animal,folder}...
            ,locImaging{animal,folder},FaceHighImaging{animal,folder},FaceLowImaging{animal,folder},locImaging_Shuffled{animal,folder},FaceHighImaging_Shuffled{animal,folder},FaceLowImaging_Shuffled{animal,folder}]...
                =stateCrossCorr2(wheelOn_final1,wheelOff_final1,Face_HighArousal_On_final1,Face_HighArousal_Off_final1,Face_LowArousal_On_final1,Face_LowArousal_Off_final1,imaging_time, dFoF_parcells,dFoF_parcells_Shuffled,params.CorrLags,names,params.signalsExtraction);
            
            % run Pearson's correlations between areas on data from three states  
            elseif strcmp(params.corrType,'pearsons')
             [loccorrMat{animal,folder},loccorrMat_Shuffled{animal,folder},FaceHighcorrMat{animal,folder},FaceHighcorrMat_Shuffled{animal,folder},FaceLowcorrMat{animal,folder},FaceLowcorrMat_Shuffled{animal,folder}...
            ,locImaging{animal,folder},FaceHighImaging{animal,folder},FaceLowImaging{animal,folder},locImaging_Shuffled{animal,folder},FaceHighImaging_Shuffled{animal,folder},FaceLowImaging_Shuffled{animal,folder}]...
                =statePearsonCorrelation(wheelOn_final1,wheelOff_final1,Face_HighArousal_On_final1,Face_HighArousal_Off_final1,Face_LowArousal_On_final1,Face_LowArousal_Off_final1,imaging_time, dFoF_parcells,dFoF_parcells_Shuffled,names,params.signalsExtraction);
            end 
        
        end       
    end
    %average correlation matrices across trials concatenated across sessions within an animal
    indivFigureFolder=fullfile(figuresFolder,MiceAnalyze{animal});
    names=fieldnames(dFoF_parcells);
    if strcmp(params.signalsExtraction,'RCaMP_AC'),names=[names;{'BG'}]; end
    for i=1:length(names)
        if  size(loccorrMat,1)==animal
            [hs,loccorrMatCatZScore.(names{i}){animal},FaceHighcorrMatCatZScore.(names{i}){animal},FaceLowcorrMatCatZScore.(names{i}){animal},loccorrMatCatZScoreShuffled.(names{i}){animal},FaceHighcorrMatCatZScoreShuffled.(names{i}){animal},FaceLowcorrMatCatZScoreShuffled.(names{i}){animal},...
            loccorrMatCat.(names{i}){animal},FaceHighcorrMatCat.(names{i}){animal},FaceLowcorrMatCat.(names{i}){animal},loccorrMatCatShuffled.(names{i}){animal},FaceHighcorrMatCatShuffled.(names{i}){animal}, FaceLowcorrMatCatShuffled.(names{i}){animal},PermLocFaceH,PermLocFaceL,PermFaceHFaceL]...
            =averageCorrMatrixSessions_withPermutation2(loccorrMat(animal,:),FaceHighcorrMat(animal,:),FaceLowcorrMat(animal,:),loccorrMat_Shuffled(animal,:),FaceHighcorrMat_Shuffled(animal,:),FaceLowcorrMat_Shuffled(animal,:),names{i},'Locomotion','FaceHigh','FaceLow',parcellnames,params.numPerms);           
        PermuteLocFaceH.(names{i}).state1corrMatCatZScorePerm{animal}= PermLocFaceH.state1corrMatCatZScorePerm;   
        PermuteLocFaceH.(names{i}).state2corrMatCatZScorePerm{animal}= PermLocFaceH.state2corrMatCatZScorePerm;  
        PermuteLocFaceH.(names{i}).state1corrMatCatZScorePermShuff{animal}= PermLocFaceH.state1corrMatCatZScorePermShuff;   
        PermuteLocFaceH.(names{i}).state2corrMatCatZScorePermShuff{animal}= PermLocFaceH.state2corrMatCatZScorePermShuff;  
        
        PermuteLocFaceL.(names{i}).state1corrMatCatZScorePerm{animal}= PermLocFaceL.state1corrMatCatZScorePerm;   
        PermuteLocFaceL.(names{i}).state2corrMatCatZScorePerm{animal}= PermLocFaceL.state3corrMatCatZScorePerm;  
        PermuteLocFaceL.(names{i}).state1corrMatCatZScorePermShuff{animal}= PermLocFaceL.state1corrMatCatZScorePermShuff;   
        PermuteLocFaceL.(names{i}).state2corrMatCatZScorePermShuff{animal}= PermLocFaceL.state3corrMatCatZScorePermShuff;  
        
        PermuteFaceHFaceL.(names{i}).state1corrMatCatZScorePerm{animal}= PermFaceHFaceL.state2corrMatCatZScorePerm;   
        PermuteFaceHFaceL.(names{i}).state2corrMatCatZScorePerm{animal}= PermFaceHFaceL.state3corrMatCatZScorePerm;  
        PermuteFaceHFaceL.(names{i}).state1corrMatCatZScorePermShuff{animal}= PermFaceHFaceL.state2corrMatCatZScorePermShuff;   
        PermuteFaceHFaceL.(names{i}).state2corrMatCatZScorePermShuff{animal}= PermFaceHFaceL.state3corrMatCatZScorePermShuff;  
        
        saveas(hs,fullfile(indivFigureFolder,strcat(names{i},'LocomotionvsHighFacevsLowFaceCorrMatrix')));
        end
    end    
    close all;
end
%save some individual mouse output
save (fullfile(figuresFolder,'IndivMouseOutput'),'FaceHighImaging','FaceLowImaging','locImaging','FaceHighcorrMat', 'FaceLowcorrMat','loccorrMat','loccorrMatCatZScore','FaceHighcorrMatCatZScore','FaceLowcorrMatCatZScore','CompData')
%% population averages, take a mean across animals, first mean of fisher's z , then backtransform z to pearson's r for display purposes
names=fieldnames(loccorrMatCatZScore);
for i=1:length(names)
        % do stats on locomotion vs face high, locomotion vs face low, and face high vs face low
        figFolder=fullfile(figuresFolder,'LocomotionFace'); if ~exist(figFolder),mkdir(figFolder), end
        [output.(names{i}).meanLoc,output.(names{i}).meanFaceH,output.(names{i}).meanFaceL,output.(names{i}).LocFaceHmeanDiffMat,output.(names{i}).LocFaceLmeanDiffMat,output.(names{i}).FaceHFaceLmeanDiffMat,...
        output.(names{i}).AllParLocFaceH, output.(names{i}).AllParLocFaceL, output.(names{i}).AllParFaceHFaceL]= statsCorrelation_Combined(loccorrMatCatZScore.(names{i}),FaceHighcorrMatCatZScore.(names{i}),...
        FaceLowcorrMatCatZScore.(names{i}),PermuteLocFaceH.(names{i}).state1corrMatCatZScorePerm,PermuteLocFaceH.(names{i}).state2corrMatCatZScorePerm,PermuteLocFaceL.(names{i}).state1corrMatCatZScorePerm,...
        PermuteLocFaceL.(names{i}).state2corrMatCatZScorePerm,PermuteFaceHFaceL.(names{i}).state1corrMatCatZScorePerm,PermuteFaceHFaceL.(names{i}).state2corrMatCatZScorePerm,...
         names{i},'loc','FaceHigh','FaceLow',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
    
    
        figFolder=fullfile(figuresFolder,'LocomotionFaceShuffled'); if ~exist(figFolder),mkdir(figFolder), end
        [output.(names{i}).meanLocShuff,output.(names{i}).meanFaceHShuff,output.(names{i}).meanFaceLShuff,output.(names{i}).LocFaceHmeanDiffMatShuff,output.(names{i}).LocFaceLmeanDiffMatShuff,...
        output.(names{i}).FaceHFaceLmeanDiffMatShuff,output.(names{i}).AllParLocFaceHShuff, output.(names{i}).AllParLocFaceLShuff, output.(names{i}).AllParFaceHFaceLShuff]= statsCorrelation_Combined(loccorrMatCatZScoreShuffled.(names{i}),...
        FaceHighcorrMatCatZScoreShuffled.(names{i}),FaceLowcorrMatCatZScoreShuffled.(names{i}),PermuteLocFaceH.(names{i}).state1corrMatCatZScorePermShuff,PermuteLocFaceH.(names{i}).state2corrMatCatZScorePermShuff,...
        PermuteLocFaceL.(names{i}).state1corrMatCatZScorePermShuff,PermuteLocFaceL.(names{i}).state2corrMatCatZScorePermShuff,PermuteFaceHFaceL.(names{i}).state1corrMatCatZScorePermShuff,...,
        PermuteFaceHFaceL.(names{i}).state2corrMatCatZScorePermShuff,names{i},'locShuff','FaceHighShuff','FaceLowShuff',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);    
end

%save output
save(fullfile(figuresFolder,'SummaryData.mat'),'output');

%% get individual correlation matrices 
names=fieldnames(loccorrMatCatZScore);
for i=1:length(names)
        % do stats on locomotion vs face high, locomotion vs face low, and face high vs face low
        figFolder=fullfile(figuresFolder,'LocomotionFace'); if ~exist(figFolder),mkdir(figFolder), end
        makeIndivCorrMatrixPlots(loccorrMatCatZScore.(names{i}),FaceHighcorrMatCatZScore.(names{i}),FaceLowcorrMatCatZScore.(names{i}),names{i},'loc','FaceHigh','FaceLow',parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figFolder);
end 

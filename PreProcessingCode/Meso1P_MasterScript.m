%% preprocess raw tif (or cxd) video files and extract uv corrected df/f 
%% Lohani and Moberly et al. 2022
addpath(genpath('.\Functions'));% where functions are located 
%% input and output 
mainDir='R:\AnimalProjects\Widefield\Data\Test\20250721\Test_20250721_Run2'; % where raw tiff files  are stored 
TmpFile='C:\Temp\Data'; %temporary folder location in local drive 
outputFolder='R:\AnimalProjects\Widefield\Preprocessed\Test\20250721\Test_20250721_Run2'; % main output folder 
 
%% user-selected input parameters 
params.fsspike2=5000;% spike2 sampling frequency 
params.fsimaging=62; % imaging sampling frequency
params.pupilSR=30; %pupil camera sampling frequncy
params.batchSize = 4081;% load tif images in batches of this number 
params.batches2load=10000;% max batches of tif to load 
params.regressUV=1; %option of regressing every pixel by UV to remove 
params.patchSize=14; % patch size to use for hemodynamic spatial regression , patch size of 14 corresponds to 9X9 pixel squared area 
params.moveLocal=1;%if you want to move the tiff file to a local drive before reading for speed
params.ImageFormat='tif';%choose between cxd or tif 

%grabs/rcamp (two cameras, RCaMP_AC) vs grabs (one camera, blueuv) only
params.signalsExtraction.firstCaFrame = 1; %1 if first frame blue 
params.signalsExtraction.sigs = 'blueuv';% 'blueuv' or 'RCaMP_AC'
params.signalsExtraction.blueuvRatio = 1;
params.resizeRatio = 1; %downsample images 
params.angle2rotate= 180; %rotation angle; 
% detrending filter params
params.deterend.filtLen = 1000;
params.deterend.filtcutoff = 0.001;
params.deterend.method = 'FIR';

%event params 
params.visStimAn=0; % 1 if vis stim are presented
params.airpuffAn=0; % 1 if airpuffs or electrical stim are presented 
params.visDur=2; %visual stimulus duration 
params.visITI=5; % ITI for visual stim trilas
params.minRunDuration=5; %minimum run duration in seconds
params.minSitDuration=10; %minimum sit duration in seconds 
params.ITITime=5; %  dead time after events in seconds  

%% choose spike2 channels for various triggers and signals 
channels.BLUE = 2;%blue led 
channels.UV = 3;%uv led 
channels.FRAMETICKS = 1;%green/uv mesocam triggers
channels.PHOTO_DIODE = [];%visual stim trigger
channels.WHEEL = 1; % wheel signal 
channels.AIR_PUFF = []; %this is either airpuff channel or electrical stim channel 
channels.PUPIL_CAMERA = 4;%pupil camera triggers
channels.RED_MESO = [];%red mesocam trigger
channels.GREEN=[];%green LED
channels.EEG = [];%eeg continous signal 

%% get smrx and visual stim files 
binfilelist = (dir(fullfile(mainDir, '*.bin')));
dataSmrxFile=fullfile(mainDir,binfilelist.name); 
%dataVisStimList=(dir(fullfile(mainDir, '*.csv')));
%dataVisStimFile=fullfile(mainDir,dataVisStimList.name); 

subStr = extractAfter(mainDir,'Data\');
outputPath=fullfile(outputFolder,subStr);

%% preprocessing step 1 (load tiffs and separate colors)
[sigsMov, R,C]= tiffExtractionPreProcessing1(mainDir, outputPath, TmpFile, params);

%% preprocessing step 2(register two cameras, detrend data, calculate dF/F)
[dFoF,dFoF_parcells, R,C]= tiffExtractionPreProcessing2(mainDir,outputPath,sigsMov,R,C,params);

%% preprocessing step 3 (hemodynamic correction with uv signal)-- recommend doing it on high perfomrance computing clusters 
if params.regressUV
 [dFoF,dFoF_parcells, R,C]= RegressionProcessing(dFoF,R,C,params,outputPath) ; 
end

%% preprprocessing step 4 (smrx file extraction) 
smrxExtractionPreProcessing(mainDir,dataSmrxFile,dataVisStimFile,outputPath,channels,dFoF,dFoF_parcells,R,C,params)

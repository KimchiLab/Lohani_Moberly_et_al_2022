%% This code adjusts the proc outputs from faceMap (with ROIs for whole face, whisker and pupil) for all sessions by flipping the face and whisker PC1 signs 
%  based on skewness sign and normalizes all data 
%% inputs and outputs 
outputFolder='W:\GRABS_Data\AnalyzedLocalPharmacology';%input and output folder 
MiceAnalyze=[{'cn1\'},{'cn2\'},{'cn3\'}]; 
%% loop through each mouse folder 
for z=1:length(MiceAnalyze)
mainDir1 =fullfile(outputFolder,MiceAnalyze{z});     
folders=dir(mainDir1); 
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
DirFolders= folders(dirFlags);
Flag=~contains({DirFolders.name},'Extra','IgnoreCase',true);
DirFolders=DirFolders(Flag);   
for i =1: length(DirFolders)   
clearvars -except outputFolder z MiceAnalyze mainDir1 folders dirFlags DirFolders Flag DirFolders i 
currfolder=fullfile(mainDir1, DirFolders(i).name);
files=dir(currfolder); 
for t=1:length(files)
if contains(files(t).name,'proc','IgnoreCase',true)
 files=files(t);break; 
end 
end 
if isempty(files), continue, end; 
%load proc data 
load(fullfile(currfolder,files.name))
pupil=proc.pupil.area; 
face=proc.motSVD{1,1}(:,1); 
whisker=proc.motSVD{1,2}(:,1); 
 
%check skewness and flip the sign if face or whisker is negatively skewed 
faceskew=(skewness(face)); 
if sign(faceskew)==-1
    face=-face; 
end 
whiskerskew=(skewness(whisker));
if sign(whiskerskew)==-1
    whisker=-whisker; 
end 

%normalize pupil, face and PC1 
Ms = sort(pupil,'ascend');% Sort asending along time dimension
F0Vals = Ms(1:ceil(length(Ms)*0.1)); % lower 10% of the values
MeanF0=mean(F0Vals);
normPupil=(pupil-MeanF0)./std(F0Vals); 
Ms = sort(face,'ascend');% Sort asending along time dimension
F0Vals = Ms(1:ceil(length(Ms)*0.1)); % lower 10% of the values
MeanF0=mean(F0Vals);
normFace=(face-MeanF0)./std(F0Vals); 
Ms = sort(whisker,'ascend');% Sort asending along time dimension
F0Vals = Ms(1:ceil(length(Ms)*0.1)); % lower 10% of the values
MeanF0=mean(F0Vals);
normWhisker=(whisker-MeanF0)./std(F0Vals); 

figure; 
subplot(3,1,1); plot(normPupil); title('NormPupil'); 
subplot(3,1,2); plot(normFace);title('NormFacePC1_Corr'); 
subplot(3,1,3); plot(normWhisker);title('NormWhiskerPC1_Corr'); 

%save the adjusted proc outputs
proc.output.pupilNorm=normPupil; 
proc.output.facePC1CorrNorm=normFace; 
proc.output.whiskerPC1CorrNorm=normWhisker;

save(fullfile(currfolder,files.name),'proc'); close all; 
end 
end 

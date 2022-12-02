function [Color_data]=makePopAverages(Color_data,Comb_Color_dataTrials, eventName,NormType,BrainReg,params)
%summarize data and do stats
Color_data.(eventName).(NormType).(BrainReg)=cell2mat(cellfun(@(x) nanmean(x,2),Comb_Color_dataTrials.(eventName).(strcat(NormType,'_',BrainReg,'_dffIntp')),'UniformOutput',false));
%take mean in moving window do stats (animals)
numAnimals=size(Color_data.(eventName).(NormType).(BrainReg),2); frameWin=params.fsimaging*params.statWin; numFrames=size(Color_data.(eventName).(NormType).(BrainReg),1)/frameWin;
maxout=zeros(numAnimals,numFrames);
meanout=zeros(numAnimals,numFrames);
for j=1:numAnimals
    x=Color_data.(eventName).(NormType).(BrainReg)(:,j);
    maxout(j,:) = accumarray(ceil((1:size(x,1))/frameWin)',x(:),[],@max);
    meanout(j,:) = accumarray(ceil((1:size(x,1))/frameWin)',x(:),[],@mean);
end
Color_data.(eventName).(NormType).(strcat(BrainReg,'_MeanWin'))=meanout;
Color_data.(eventName).(NormType).(strcat(BrainReg,'_MaxWin'))=maxout;
%get trials for one example mouse only
Color_data.(eventName).(NormType).(strcat(BrainReg,'_ExampleMouse'))=Comb_Color_dataTrials.(eventName).(strcat(NormType,'_',BrainReg,'_dffIntp')){params.Mouse};
end
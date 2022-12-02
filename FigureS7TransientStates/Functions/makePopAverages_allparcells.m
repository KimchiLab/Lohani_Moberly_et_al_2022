function [figure1,figure2, Color_data]=makePopAverages_allparcells(Color_data,Comb_Color_dataTrials, eventName,NormType,params,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,parcellnames,ylimBrainMap)
%summarize data and do stats
BrainReg='allParcells';
Color_data.(eventName).(NormType).(BrainReg)=cell2mat(cellfun(@(x) nanmean(x,2),Comb_Color_dataTrials.(eventName).(strcat(NormType,'_parcellsDFF')),'UniformOutput',false));
%take mean in moving window do stats (animals)
numAnimals=size(Color_data.(eventName).(NormType).(BrainReg),2); frameWin=params.fsimaging*params.statWin; numFrames=size(Color_data.(eventName).(NormType).(BrainReg),1)/frameWin;
numParcells=size(Color_data.(eventName).(NormType).(BrainReg),3);
maxout=zeros(numAnimals,numFrames,numParcells);
meanout=zeros(numAnimals,numFrames,numParcells);
for j=1:numAnimals
    for k=1:numParcells
        x=Color_data.(eventName).(NormType).(BrainReg)(:,j,k);
        maxout(j,:,k) = accumarray(ceil((1:size(x,1))/frameWin)',x(:),[],@max);
        meanout(j,:,k) = accumarray(ceil((1:size(x,1))/frameWin)',x(:),[],@mean);
    end
end

%get population average across animals in the window [0 1] following transition point
frame=((params.preEventWin+params.statTime)*params.statWin+1);
popAve_meanWin=squeeze(nanmean(meanout(:,frame,:),1));
popAve_maxWin=squeeze(nanmean(maxout(:,frame,:),1));
alldata_meanWin=squeeze(meanout(:,frame,:));
alldata_maxWin=squeeze(maxout(:,frame,:));
%plot max win data (population average) as whole brain maps
[figure1]=plot_parcel_correlation_brainmap(popAve_maxWin(CombinedParcellIdx),'parula',ylimBrainMap,CombinedParcellIdx,V1Idx,M2Idx,S1bIdx); %brain color map
%plot animal by animal max win data as bar graph only for left hemisphere
figure2=figure; numParcells=length(CombinedParcellIdx)/2;
data=alldata_maxWin(:,CombinedParcellIdx(2:2:end))';
bar(popAve_maxWin(CombinedParcellIdx(2:2:end)),'FaceColor','none');hold on; scatter(repmat(1:numParcells,1,numAnimals),data(:),10,'k','filled');ylabel('Peak Z');
xticks(1:numParcells); xticklabels(parcellnames(2:2:end));xtickangle(90); ylim([-0.1 1]);
box off; set(gca,'TickDir','out');

Color_data.(eventName).(NormType).(strcat(BrainReg,'_MeanWin'))=alldata_meanWin;
Color_data.(eventName).(NormType).(strcat(BrainReg,'_MaxWin'))=alldata_maxWin;

end


function makeIndivCorrMatrixPlots(dataMatState1,dataMatState2,dataMatState3,fieldname,BehavStat1,BehavStat2,BehavStat3,parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figuresFolder)
% This function make individual correlation matrix plots
%inputs: dataMatState1 for state 1 as Zscore corr matrix for each animal, will be converted to r-value for display
%dataMatState2 for state 2 as Zscore corr matrix for each animal, will be converted to r-value for display
%dataMatState3 for state 3 as Zscore corr matrix for each animal, will be converted to r-value for display
%fieldname: blue, uv, green or BG
%BehavStat1:name of State1
%BehavStat2: name of State2
%BehavStat3: name of State3
%parcellnames: names of parcells
%CombinedParcellIdx: indices of parcells we are using
%V1Ids,S1bIdx,M2Idx: parcell idx for the three specific parcells
%figuresFolder: output folder for figures and data

%%
dataMat1=cat(3,dataMatState1{:});
dataMat2=cat(3,dataMatState2{:});
dataMat3=cat(3,dataMatState3{:});

diffMat12=dataMat1-dataMat2;
diffMat13=dataMat1-dataMat3;
diffMat23=dataMat2-dataMat3;

numparcells=size(dataMat1,1);
numanimals=size(dataMat1,3);
if ~strcmp(fieldname,'BG')       %if the field is not blue green, get the lower triangular part of matrix and make correlation matrices
    figure1=figure;
    startPlot=0;
    for animal=1:numanimals
        meanState1=ifisherz(dataMat1(:,:,animal));
        meanState2=ifisherz(dataMat2(:,:,animal));
        meanState3=ifisherz(dataMat3(:,:,animal));
        
        meanDiffMat12=ifisherz(diffMat12(:,:,animal));
        meanDiffMat13=ifisherz(diffMat13(:,:,animal));
        meanDiffMat23=ifisherz(diffMat23(:,:,animal));
        %set the upper triangular part of matrix to zero
        tmp=ones(numparcells,numparcells);
        tmp=tril(tmp-diag(diag(tmp)));
        [Idx]=find(tmp==0);
        
        meanState1(Idx)=nan;
        meanState2(Idx)=nan;
        meanState3(Idx)=nan;
        meanDiffMat12(Idx)=nan;
        meanDiffMat13(Idx)=nan;
        meanDiffMat23(Idx)=nan;
        
        %make mean correlation matrix plots
        ax(1+startPlot)=subplot(numanimals,6,1+startPlot);
        imagesc(meanState1,'AlphaData',~isnan(meanState1)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); line(xlim,ylim,'Color',[0 0 0]); %regular correlation matrix
        caxis([0 1]); title(strcat(fieldname,BehavStat1,'Animal#',num2str(animal))); box off;
        xtickangle(90); colormap(ax(1+startPlot),parula); set(gca,'TickDir','out');
        
        xticklabels(parcellnames); yticklabels(parcellnames(:)); colorbar;
        ax(2+startPlot)=subplot(numanimals,6,2+startPlot);
        imagesc(meanState2,'AlphaData',~isnan(meanState2));xticks(1:length(parcellnames));yticks(1:length(parcellnames)); line(xlim,ylim,'Color',[0 0 0]);
        caxis([0 1]); title(strcat(fieldname,BehavStat2,'Animal#',num2str(animal))); box off;
        xtickangle(90); colormap(ax(2+startPlot),parula); set(gca,'TickDir','out');
        xticklabels(parcellnames); yticklabels(parcellnames(:)); colorbar;
        
        ax(3+startPlot)=subplot(numanimals,6,3+startPlot);
        imagesc(meanState3,'AlphaData',~isnan(meanState3));xticks(1:length(parcellnames));yticks(1:length(parcellnames)); line(xlim,ylim,'Color',[0 0 0]);
        caxis([0 1]); title(strcat(fieldname,BehavStat3,'Animal#',num2str(animal))); box off;
        xtickangle(90); colormap(ax(3+startPlot),parula); set(gca,'TickDir','out');
        xticklabels(parcellnames); yticklabels(parcellnames(:)); colorbar;
        
        %make difference correlation matrix plot
        ax(4+startPlot)=subplot(numanimals,6,4+startPlot);
        imagesc(meanDiffMat12,'AlphaData',~isnan(meanDiffMat12)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); line(xlim,ylim,'Color',[0 0 0]);
        caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat1,'-',BehavStat2,'Diff-Animal#',num2str(animal))); box off;
        xtickangle(90); colormap(ax(4+startPlot),redblue);set(gca,'TickDir','out');
        xticklabels(parcellnames); yticklabels(parcellnames(:)); colorbar;
        
        ax(5+startPlot)=subplot(numanimals,6,5+startPlot);
        imagesc(meanDiffMat13,'AlphaData',~isnan(meanDiffMat13)); xticks(1:length(parcellnames));yticks(1:length(parcellnames));line(xlim,ylim,'Color',[0 0 0]);
        caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat1,'-',BehavStat3,'Diff-Animal#',num2str(animal))); box off;
        xtickangle(90); colormap(ax(5+startPlot),redblue); set(gca,'TickDir','out');
        xticklabels(parcellnames); yticklabels(parcellnames(:)); colorbar;
        
        ax(6+startPlot)=subplot(numanimals,6,6+startPlot);
        imagesc(meanDiffMat23,'AlphaData',~isnan(meanDiffMat23)); xticks(1:length(parcellnames));yticks(1:length(parcellnames));line(xlim,ylim,'Color',[0 0 0]);
        caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat2,'-',BehavStat3,'Diff-Animal#',num2str(animal))); box off;
        xtickangle(90); colormap(ax(6+startPlot),redblue);set(gca,'TickDir','out');
        xticklabels(parcellnames); yticklabels(parcellnames(:)); colorbar;
        
        startPlot=startPlot+6;
    end
    %save figures
    saveas(figure1, fullfile(figuresFolder,strcat(fieldname,'CorrelationsIndivAnimalState')));
    saveas(figure1, fullfile(figuresFolder,strcat(fieldname,'CorrelationsIndivAnimalState.pdf')));
    
else %if the field is blue green, only focus on diagonal elements of the matrix and plot correlations as brain color maps
    for animal=1:numanimals
        meanState1=ifisherz(dataMat1(:,:,animal));
        meanState2=ifisherz(dataMat2(:,:,animal));
        meanState3=ifisherz(dataMat3(:,:,animal));
        meanDiffMat12=ifisherz(diffMat12(:,:,animal));
        meanDiffMat23=ifisherz(diffMat23(:,:,animal));
        
        dataState1=diag(meanState1);
        dataState2=diag(meanState2);
        dataState3=diag(meanState3);
        dataStateDiff12=diag(meanDiffMat12);
        dataStateDiff23=diag(meanDiffMat23);
    
    [figure1]=plot_parcel_correlation_brainmap(dataState1,'parula',[0 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    set(gca,'XDir','reverse');title(strcat(fieldname,BehavStat1,'Animal#',num2str(animal)));
    
    
    [figure2]=plot_parcel_correlation_brainmap(dataState2,'parula',[0 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    set(gca,'XDir','reverse');title(strcat(fieldname,BehavStat2,'Animal#',num2str(animal)));
    
  
    [figure3]=plot_parcel_correlation_brainmap(dataState3,'parula',[0 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    set(gca,'XDir','reverse');title(strcat(fieldname,BehavStat3,'Animal#',num2str(animal)));
    
    
    [figure4]=plot_parcel_correlation_brainmap(dataStateDiff12,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    set(gca,'XDir','reverse');title(strcat(fieldname,BehavStat1,'-',BehavStat2,'Diff-Animal#',num2str(animal)));
    
    
    [figure5]=plot_parcel_correlation_brainmap(dataStateDiff23,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    set(gca,'XDir','reverse');title(strcat(fieldname,BehavStat2,'-',BehavStat3,'Diff-Animal#',num2str(animal)));
    
    saveas(figure1, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'Animal#',num2str(animal))));
    saveas(figure1, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'Animal#',num2str(animal),'.pdf')));
    saveas(figure2, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'Animal#',num2str(animal))));
    saveas(figure2, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'Animal#',num2str(animal),'.pdf')));
    saveas(figure3, fullfile(figuresFolder,strcat(fieldname,BehavStat3,'Animal#',num2str(animal))));
    saveas(figure3, fullfile(figuresFolder,strcat(fieldname,BehavStat3,'Animal#',num2str(animal),'.pdf')));
    saveas(figure4, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat2,'Diff-Animal#',num2str(animal))));
    saveas(figure4, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat2,'Diff-Animal#',num2str(animal),'.pdf')));
    saveas(figure5, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-',BehavStat3,'Diff-Animal#',num2str(animal))));
    saveas(figure5, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-',BehavStat3,'Diff-Animal#',num2str(animal),'.pdf')));
    end 
end
end


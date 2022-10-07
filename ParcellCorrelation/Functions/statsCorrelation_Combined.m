function [meanState1,meanState2,meanState3,meanDiffMat12,meanDiffMat13,meanDiffMat23,AllPar12,AllPar13,AllPar23]= statsCorrelation_Combined(dataMatState1,dataMatState2,dataMatState3, dataMat12State1Perm,dataMat12State2Perm,...
    dataMat13State1Perm,dataMat13State2Perm,dataMat23State1Perm,dataMat23State2Perm,fieldname,BehavStat1,BehavStat2,BehavStat3,parcellnames,CombinedParcellIdx,V1Idx,S1bIdx,M2Idx,figuresFolder)
% This function does stats using paired permutations
%inputs: dataMatState1 for state 1 as Zscore corr matrix for each animal, will be converted to r-value for display
%dataMatState2 for state 2 as Zscore corr matrix for each animal, will be converted to r-value for display
%dataMatState3 for state 3 as Zscore corr matrix for each animal, will be converted to r-value for display
%dataMat12State1Perm/dataMat12State2Perm: randomly permuted corr matrix for each animal for comparison between states 1 and 2 
%dataMat13State1Perm/dataMat13State2Perm: randomly permuted corr matrix for each animal for comparison between states 1 and 3
%dataMat23State1Perm/dataMat23State2Perm: randomly permuted corr matrix for each animal for comparison between states 2 and 3
%fieldname: blue, uv, green or BG
%BehavStat1:name of State1
%BehavStat2: name of State2
%BehavStat3: name of State3
%parcellnames: names of parcells
%CombinedParcellIdx: indices of parcells we are using
%V1Ids,S1bIdx,M2Idx: parcell idx for the three specific parcells
%figuresFolder: output folder for figures and data

%%ouputs
%meanState1/meanState2/meanState3:correlation matrix averaged across animals for State 1, State2 and State3 
%meanDiffMat12: difference correlation matrix between State1 and 2 averaged across animals
%meanDiffMat13: difference correlation matrix between State1 and 3 averaged across animals
%meanDiffMat23: difference correlation matrix between State2 and 3 averaged across animals
%AllPar12/AllPar13/AllPar23: output for all parcells, including stats for respective state comparisons(Between states 1 and 2, or states 1 and 3 or states 2 and 3) 
%%
%concatenate data across animals
dataMat1=cat(3,dataMatState1{:});
dataMat2=cat(3,dataMatState2{:});
dataMat3=cat(3,dataMatState3{:});

numparcells=size(dataMat1,1);
numanimals=size(dataMat1,3);

%make a matrix of mean of State1 and State2 by first taking a mean of Z-scores across animals and converting to r for display purpose
%also make a matrix of difference between State 1 and State2, take a mean of differences  and convert to r for display purpose
%also do pairwise stats between State1 and State2  on each cell and show the difference matrix as values for significant cells only
meanState1=ifisherz(nanmean(dataMat1,3));
meanState2=ifisherz(nanmean(dataMat2,3));
meanState3=ifisherz(nanmean(dataMat3,3));
diffMat12=dataMat1-dataMat2;
meanDiffMat12=ifisherz(nanmean(diffMat12,3));

diffMat13=dataMat1-dataMat3;
meanDiffMat13=ifisherz(nanmean(diffMat13,3));

diffMat23=dataMat2-dataMat3;
meanDiffMat23=ifisherz(nanmean(diffMat23,3));

if ~strcmp(fieldname,'BG')       %if the field is not blue green, get the lower triangular part of matrix and make correlation matrices
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
    figure1=figure;
    cmap=parula; 
    imagesc(meanState1,'AlphaData',~isnan(meanState1)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:));  line(xlim,ylim,'Color',[0 0 0]); %regular correlation matrix
    caxis([0 1]); title(strcat(fieldname,BehavStat1,'-Mean')); box off;
    xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');
    
    figure2=figure;
    imagesc(meanState2,'AlphaData',~isnan(meanState2));xticks(1:length(parcellnames));yticks(1:length(parcellnames));xticklabels(parcellnames); yticklabels(parcellnames(:)); line(xlim,ylim,'Color',[0 0 0]);
    caxis([0 1]); title(strcat(fieldname,BehavStat2,'-Mean')); box off;
    xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');
    
    figure3=figure;
    imagesc(meanState3,'AlphaData',~isnan(meanState3));xticks(1:length(parcellnames));yticks(1:length(parcellnames));xticklabels(parcellnames); yticklabels(parcellnames(:)); line(xlim,ylim,'Color',[0 0 0]);
    caxis([0 1]); title(strcat(fieldname,BehavStat3,'-Mean')); box off;
    xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');
    
    %make difference correlation matrix plot
    figure4=figure;
    imagesc(meanDiffMat12,'AlphaData',~isnan(meanDiffMat12)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:));line(xlim,ylim,'Color',[0 0 0]);
    caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiff')); box off;
    xtickangle(90); colormap(redblue); colorbar;set(gca,'TickDir','out');
    
    figure5=figure;
    imagesc(meanDiffMat13,'AlphaData',~isnan(meanDiffMat13)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:));line(xlim,ylim,'Color',[0 0 0]);
    caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiff')); box off;
    xtickangle(90); colormap(redblue); colorbar;set(gca,'TickDir','out');
    
    figure6=figure;
    imagesc(meanDiffMat23,'AlphaData',~isnan(meanDiffMat23)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:));line(xlim,ylim,'Color',[0 0 0]);
    caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiff')); box off;
    xtickangle(90); colormap(redblue); colorbar;set(gca,'TickDir','out');
    %% do paired permutation test followed by multiple comparison correction using data previoulsy randomized/permuted by shuffling state epoch labels within each animal
    [AllPar12]=permutePval(numparcells,numanimals,dataMat12State1Perm,dataMat12State2Perm,dataMat1,dataMat2);
    [AllPar13]=permutePval(numparcells,numanimals,dataMat13State1Perm,dataMat13State2Perm,dataMat1,dataMat3);
    [AllPar23]=permutePval(numparcells,numanimals,dataMat23State1Perm,dataMat23State2Perm,dataMat2,dataMat3);

    %show difference correlation matrix with values for significant cells only, non significant ones are set to zero
    tmp=ones(numparcells,numparcells);
    tmp=tril(tmp-diag(diag(tmp)));
    [Idx]=find(tmp==1);  
    
    figure7=figure;
    meanDiffMat1=meanDiffMat12;
    indices=find(AllPar12.IndivCell.permBHC.hVal==0);%find non-significant cells
    indices=intersect(indices,Idx); %only get cells in the lower triangle 
    meanDiffMat1(indices)=0;%set non-significant cells to zero
    cmap=[redblue(512)]; cmap(257,:)=[211/255 211/255 211/255]; 
    imagesc(meanDiffMat1,'AlphaData',~isnan(meanDiffMat1)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:));line(xlim,ylim,'Color',[0 0 0]);
    caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiff-OnlySigParcells')); box off;
    xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');
    
    figure8=figure;
    meanDiffMat2=meanDiffMat13;
    indices=find(AllPar13.IndivCell.permBHC.hVal==0);%find non-significant cells
    indices=intersect(indices,Idx); %only get cells in the lower triangle 
    meanDiffMat2(indices)=0;%set non-significant cells to zero
    imagesc(meanDiffMat2,'AlphaData',~isnan(meanDiffMat2)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:));line(xlim,ylim,'Color',[0 0 0]);
    caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiff-OnlySigParcells')); box off;
    xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');
    
    figure9=figure;
    meanDiffMat3=meanDiffMat23;
    indices=find(AllPar23.IndivCell.permBHC.hVal==0);%find non-significant cells
    indices=intersect(indices,Idx); %only get cells in the lower triangle 
    meanDiffMat3(indices)=0;%set non-significant cells to zero
    imagesc(meanDiffMat3,'AlphaData',~isnan(meanDiffMat3)); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:));line(xlim,ylim,'Color',[0 0 0]);
    caxis([-0.5 0.5]); title(strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiff-OnlySigParcells')); box off;
    xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');
    
else %if the field is blue green, only focus on diagonal elements of the matrix and plot correlations as brain color maps
    dataState1=diag(meanState1);
    dataState2=diag(meanState2);
    dataState3=diag(meanState3);
    dataStateDiff12=diag(meanDiffMat12);
    dataStateDiff13=diag(meanDiffMat13);
    dataStateDiff23=diag(meanDiffMat23);
      
    [figure1]=plot_parcel_correlation_brainmap(dataState1,'parula',[0 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    title(strcat(fieldname,BehavStat1,'-Mean'));
    [figure2]=plot_parcel_correlation_brainmap(dataState2,'parula',[0 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    title(strcat(fieldname,BehavStat2,'-Mean'));
    [figure3]=plot_parcel_correlation_brainmap(dataState3,'parula',[0 0.7],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    title(strcat(fieldname,BehavStat3,'-Mean'));
    [figure4]=plot_parcel_correlation_brainmap(dataStateDiff12,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    title(strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiff'));
    [figure5]=plot_parcel_correlation_brainmap(dataStateDiff13,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    title(strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiff'));
    [figure6]=plot_parcel_correlation_brainmap(dataStateDiff23,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx);
    title(strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiff'));
  
    %% do paired permutation test followed by multiple comparison correction using data previoulsy randomized/permuted by shuffline state epoch labels within each animal  
    [AllPar12]=permutePval_BG(numparcells,numanimals,dataMat12State1Perm,dataMat12State2Perm,dataMat1,dataMat2);
    [AllPar13]=permutePval_BG(numparcells,numanimals,dataMat13State1Perm,dataMat13State2Perm,dataMat1,dataMat3);
    [AllPar23]=permutePval_BG(numparcells,numanimals,dataMat23State1Perm,dataMat23State2Perm,dataMat2,dataMat3);
         
    %show difference correlation matrix with values for significant cells only, non significant ones are set to zero
    dataStateDiff1=dataStateDiff12;
    [figure7]=plot_parcel_correlation_brainmap(dataStateDiff1,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx,AllPar12.IndivCell.permBHC.PVal);
    title(strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiff-sigparcellsOnly'));
    
    dataStateDiff2=dataStateDiff13;
    [figure8]=plot_parcel_correlation_brainmap(dataStateDiff2,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx,AllPar13.IndivCell.permBHC.PVal);
    title(strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiff-sigparcellsOnly'));
    
    dataStateDiff3=dataStateDiff23;
    [figure9]=plot_parcel_correlation_brainmap(dataStateDiff3,'redblue',[-0.6 0.6],CombinedParcellIdx,V1Idx,M2Idx,S1bIdx,AllPar23.IndivCell.permBHC.PVal);
    title(strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiff-sigparcellsOnly'));    
end
%save figures
saveas(figure1, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-Mean')));
saveas(figure1, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-Mean.pdf')));
saveas(figure2, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-Mean')));
saveas(figure2, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-Mean.pdf')));
saveas(figure3, fullfile(figuresFolder,strcat(fieldname,BehavStat3,'-Mean')));
saveas(figure3, fullfile(figuresFolder,strcat(fieldname,BehavStat3,'-Mean.pdf')));

saveas(figure4, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiff')));
saveas(figure4, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiff.pdf')));
saveas(figure5, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiff')));
saveas(figure5, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiff.pdf')));
saveas(figure6, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiff')));
saveas(figure6, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiff.pdf')));

saveas(figure7, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiffSigParcellsOnly')));
saveas(figure7, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat2,'-MeanDiffSigParcellsOnly.pdf')));
saveas(figure8, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiffSigParcellsOnly')));
saveas(figure8, fullfile(figuresFolder,strcat(fieldname,BehavStat1,'-',BehavStat3,'-MeanDiffSigParcellsOnly.pdf')));
saveas(figure9, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiffSigParcellsOnly')));
saveas(figure9, fullfile(figuresFolder,strcat(fieldname,BehavStat2,'-',BehavStat3,'-MeanDiffSigParcellsOnly.pdf')));

end

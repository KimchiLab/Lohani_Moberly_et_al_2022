function [hs,state1corrMatCatZScore,state2corrMatCatZScore,state3corrMatCatZScore,state1corrMatCatZScoreShuffled,state2corrMatCatZScoreShuffled,state3corrMatCatZScoreShuffled,...
    state1corrMatCat,state2corrMatCat,state3corrMatCat,state1corrMatCatShuffled,state2corrMatCatShuffled,state3corrMatCatShuffled,PermuteState12,PermuteState13,PermuteState23]...
    =averageCorrMatrixSessions_withPermutation2(currDataState1,currDataState2,currDataState3,currDataState1Shuffled,currDataState2Shuffled,currDataState3Shuffled,name,State1Name,State2Name,State3Name,parcellnames,numPerms)
%% This function concatenates data across sessions within an animal, performs fisher's Z-transformation on correlation matricces, takes a mean across all trials, and backtransforms to r-value for plotting
%inputs. Also performs permutation of state labels across trials and takes a mean of permuted Z-score data across all trials 
% currDataState1/currDataState2/currDataState3: correlation matrices X trials for each  session for states 1 and 2
% currDataState1Shuffled/currDataState2Shuffled/currDataState3Shuffled: correlation matrices X trials for each  session for states 1 and 2 (for shuffled data)
%name: current fieldname 
%State1Name/State2Name/State3Name: names of current states 
%parcellnames: names of parcells
%numPerms: number of permutations 
%outputs 
%hs: ouptupt figure 
%state1corrMatCatZScore,state2corrMatCatZScore,state3corrMatCatZScore: averaged and fisher's Z-transformed correlation matrix 
%state1ccorrMatCatZScoreShuffled,state2corrMatCatZScoreShuffled,state3corrMatCatZScoreShuffled: averaged and fisher's Z-transformed correlation matrix for shuffled data 
%state1corrMatCat,state2corrMatCat,state3corrMatCat: averaged and backtransformed r-value correlation matrix 
%state1corrMatCatShuffled, state2corrMatCatShuffled,state3corrMatCatShuffled:  averaged and backtransformed r-value correlation matrix for shuffled data 
%PermuteState12/PermuteState13/PermuteState23: contains permuted fisher's Z-transformed correlation matrices for raw and shuffled data (each structure contains data from permutation of epoch labels across two indicated states)  
%% written by SL, 2021
%concatenate data across sessions
newData=currDataState1(~cellfun('isempty',currDataState1));
tmp=(cellfun(@(x) x.(name), newData,'UniformOutput',false));
finalData_state1=cat(3,tmp{:});

newData=currDataState2(~cellfun('isempty',currDataState2));
tmp=(cellfun(@(x) x.(name), newData,'UniformOutput',false));
finalData_state2=cat(3,tmp{:});

newData=currDataState3(~cellfun('isempty',currDataState3));
tmp=(cellfun(@(x) x.(name), newData,'UniformOutput',false));
finalData_state3=cat(3,tmp{:});

newData=currDataState1Shuffled(~cellfun('isempty',currDataState1Shuffled));
tmp=(cellfun(@(x) x.(name), newData,'UniformOutput',false));
finalData_state1Shuff=cat(3,tmp{:});

newData=currDataState2Shuffled(~cellfun('isempty',currDataState2Shuffled));
tmp=(cellfun(@(x) x.(name), newData,'UniformOutput',false));
finalData_state2Shuff=cat(3,tmp{:}); 

newData=currDataState3Shuffled(~cellfun('isempty',currDataState3Shuffled));
tmp=(cellfun(@(x) x.(name), newData,'UniformOutput',false));
finalData_state3Shuff=cat(3,tmp{:});

%first convert correlation coefficient to Z-score (fisher's z) then take a mean
state1corrMatCatZScore=squeeze(nanmean(atanh(finalData_state1),3));
state2corrMatCatZScore=squeeze(nanmean(atanh(finalData_state2),3));
state3corrMatCatZScore=squeeze(nanmean(atanh(finalData_state3),3));
state1corrMatCatZScoreShuffled=squeeze(nanmean(atanh(finalData_state1Shuff),3));
state2corrMatCatZScoreShuffled=squeeze(nanmean(atanh(finalData_state2Shuff),3));
state3corrMatCatZScoreShuffled=squeeze(nanmean(atanh(finalData_state3Shuff),3));

%backtransform to corrleation coefficient for display purpose
state1corrMatCat=ifisherz(state1corrMatCatZScore);
state2corrMatCat=ifisherz(state2corrMatCatZScore);
state3corrMatCat=ifisherz(state3corrMatCatZScore);
state1corrMatCatShuffled=ifisherz(state1corrMatCatZScoreShuffled);
state2corrMatCatShuffled=ifisherz(state2corrMatCatZScoreShuffled);
state3corrMatCatShuffled=ifisherz(state3corrMatCatZScoreShuffled);

%plot the upper triangular part of correlation matrices
hs=figure;
L = tril(state1corrMatCat);
subplot(3,3,1); imagesc(L);caxis([0.2 1]); colormap(parula);title(strcat(name,State1Name)); box off;
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state2corrMatCat);
subplot(3,3,2);imagesc(L);caxis([0.2 1]);colormap(parula); title(strcat(name,State2Name));box off
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state3corrMatCat);
subplot(3,3,3);imagesc(L);caxis([0.2 1]);colormap(parula); title(strcat(name,State3Name));box off
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state1corrMatCatShuffled);
subplot(3,3,4); imagesc(L);caxis([0.2 1]);colormap(parula); title(strcat(name,State1Name,'Shuffled')); box off;
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state2corrMatCatShuffled);
subplot(3,3,5);imagesc(L);caxis([0.2 1]);colormap(parula); title(strcat(name,State2Name,'Shuffled'));box off
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state3corrMatCatShuffled);
subplot(3,3,6);imagesc(L);caxis([0.2 1]);colormap(parula); title(strcat(name,State3Name,'Shuffled'));box off
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state1corrMatCat-state2corrMatCat); %get the difference between state1 and state 2 
subplot(3,3,7);imagesc(L);caxis([-0.5 0.5]); colormap(redblue);title(strcat(name,State1Name,'-',State2Name));box off;
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state1corrMatCat-state3corrMatCat); %get the difference between state1 and state 3
subplot(3,3,8);imagesc(L);caxis([-0.5 0.5]); colormap(redblue);title(strcat(name,State1Name,'-',State3Name));box off;
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

L = tril(state2corrMatCat-state3corrMatCat); %get the difference between state2 and state 3
subplot(3,3,9);imagesc(L);caxis([-0.5 0.5]); colormap(redblue);title(strcat(name,State2Name,'-',State3Name));box off;
xticks(1:length(parcellnames)); xticklabels(parcellnames); yticks(1:length(parcellnames)); yticklabels(parcellnames); xtickangle(90);

%% do permutation of state labels across epochs (for comparison between state1 and state2) 
[PermuteState12.state1corrMatCatZScorePerm,PermuteState12.state2corrMatCatZScorePerm,PermuteState12.state1corrMatCatZScorePermShuff,PermuteState12.state2corrMatCatZScorePermShuff]=...
    permuteEpochLabels(finalData_state1, finalData_state2,finalData_state1Shuff,finalData_state2Shuff,numPerms);
%do permutation of state labels across epochs (for comparison between state1 and state3) 
[PermuteState13.state1corrMatCatZScorePerm,PermuteState13.state3corrMatCatZScorePerm,PermuteState13.state1corrMatCatZScorePermShuff,PermuteState13.state3corrMatCatZScorePermShuff]=...
    permuteEpochLabels(finalData_state1, finalData_state3,finalData_state1Shuff,finalData_state3Shuff,numPerms);
%do permutation of state labels across epochs (for comparison between state2 and state3) 
[PermuteState23.state2corrMatCatZScorePerm,PermuteState23.state3corrMatCatZScorePerm,PermuteState23.state2corrMatCatZScorePermShuff,PermuteState23.state3corrMatCatZScorePermShuff]=...
    permuteEpochLabels(finalData_state2, finalData_state3,finalData_state2Shuff,finalData_state3Shuff,numPerms);

end
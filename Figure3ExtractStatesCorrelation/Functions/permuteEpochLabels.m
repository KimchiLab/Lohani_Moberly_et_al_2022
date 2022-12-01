function[state1corrMatCatZScorePerm,state2corrMatCatZScorePerm,state1corrMatCatZScorePermShuff,state2corrMatCatZScorePermShuff]=...
    permuteEpochLabels(finalData_state1, finalData_state2,finalData_state1Shuff,finalData_state2Shuff,numPerms)
% do permutations of state epoch labels 
%SL,2021
nlabels=2*size(finalData_state1,3); %total number of epochs 
combinedStateData=cat(3,finalData_state1,finalData_state2); %combined raw data 
combinedStateDataShuff=cat(3,finalData_state1Shuff,finalData_state2Shuff); %combined shuffleddata 

state1corrMatCatZScorePerm=zeros(size(finalData_state1,1),size(finalData_state1,2),numPerms); 
state2corrMatCatZScorePerm=zeros(size(finalData_state2,1),size(finalData_state2,2),numPerms); 

state1corrMatCatZScorePermShuff=zeros(size(finalData_state1,1),size(finalData_state1,2),numPerms); 
state2corrMatCatZScorePermShuff=zeros(size(finalData_state2,1),size(finalData_state2,2),numPerms); 
for nPerm=1:numPerms 
    %on raw data 
    currLabels=randperm((nlabels)); %randomize epoch labels 
    finalDataPerm_state1=combinedStateData(:,:,currLabels(1:nlabels/2)); 
    finalDataPerm_state2=combinedStateData(:,:,currLabels(((nlabels/2)+1):nlabels)); 
    state1corrMatCatZScorePerm(:,:,nPerm)=squeeze(nanmean(atanh(finalDataPerm_state1),3)); %convert correlations to fisher's z and averge across randomly permuted epochs before allocating to permutation matrix 
    state2corrMatCatZScorePerm(:,:,nPerm)=squeeze(nanmean(atanh(finalDataPerm_state2),3)); %convert correlations to fisher's z and average across randomly permuted epochs before allocating to permutation matrix 
    
    %on shuffled data 
    finalDataPerm_state1Shuff=combinedStateDataShuff(:,:,currLabels(1:nlabels/2)); 
    finalDataPerm_state2Shuff=combinedStateDataShuff(:,:,currLabels(((nlabels/2)+1):nlabels)); 
    state1corrMatCatZScorePermShuff(:,:,nPerm)=squeeze(nanmean(atanh(finalDataPerm_state1Shuff),3)); %convert correlations to fisher's z and averge across randomly permuted epochs before allocating to permutation matrix 
    state2corrMatCatZScorePermShuff(:,:,nPerm)=squeeze(nanmean(atanh(finalDataPerm_state2Shuff),3)); %convert correlations to fisher's z and average across randomly permuted epochs before allocating to permutation matrix 
end
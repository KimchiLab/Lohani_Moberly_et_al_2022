function [AllPar]=permutePval_BG(numparcells,numanimals,dataMatState1Perm,dataMatState2Perm,dataMat1,dataMat2)
%This function computes p-value for comparison of state1 and state2 using permuted distribution as null distribution
%paired permutation test followed by multiple comparison correction using data previoulsy randomized/permuted by shuffling state epoch labels within each animal   
%get mean difference from permuted data 
    numPerms=size(dataMatState1Perm{1},3); %number of permutations
    meanDiffPerm=zeros(numPerms,numparcells);
    for numP=1:numPerms
        diagMatPerm=zeros(numanimals,numparcells);
        for numA=1:numanimals
            diffMatPerm=squeeze(dataMatState1Perm{numA}(:,:,numP)-dataMatState2Perm{numA}(:,:,numP));%subtract two state permuted  matrices
            diagMatPerm(numA,:)=diag(diffMatPerm);%get diagonal elements only 
        end
        meanDiffPerm(numP,:)=nanmean(diagMatPerm,1);
    end
    
    %observed mean difference 
    meanDiffData=diag(nanmean(dataMat1-dataMat2,3));
    
    %calculate p-values 
    pval=nan(1,length(meanDiffData));
    for comb=1:numparcells
        pval(comb) =  mean(abs(meanDiffPerm(:,comb))>=abs(meanDiffData(comb)));
        if pval(comb)==0
            pval(comb)=1/numPerms; %so that p-value isn't exactly zero 
        end 
    end
    hval=pval<0.05;
    AllPar.IndivCell.permUC.PVal=pval;
    AllPar.IndivCell.permUC.hVal=hval; 
    
    [hval,~,~,pval]=fdr_bh(pval,0.01,'dep'); %benjamini hotchberg correction at q =0.01 fdr 
    AllPar.IndivCell.permBHC.PVal=pval;
    AllPar.IndivCell.permBHC.hVal=hval; 
   
end 
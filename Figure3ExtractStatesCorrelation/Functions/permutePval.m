function [AllPar]=permutePval(numparcells,numanimals,dataMatState1Perm,dataMatState2Perm,dataMat1,dataMat2)
%This function computes p-value for comparison of state1 and state2 using permuted distribution as null distribution
%paired permutation test followed by multiple comparison correction using data previoulsy randomized/permuted by shuffling state epoch labels within each animal
tmp=ones(numparcells,numparcells);
tmp=tril(tmp-diag(diag(tmp)));
[Idx]=find(tmp==1);
numComp=length(Idx);

%get mean difference from permuted data
numPerms=size(dataMatState1Perm{1},3); %number of permutations
meanDiffPerm=zeros(numPerms,numComp); %variable where permuted stats are stored
for numP=1:numPerms
    MatPerm=zeros(numanimals,numComp);
    for numA=1:numanimals
        diffMatPerm=squeeze(dataMatState1Perm{numA}(:,:,numP)-dataMatState2Perm{numA}(:,:,numP));%subtract two state permuted  matrices
        MatPerm(numA,:)=diffMatPerm(Idx);%transorfm matrix to linear vector
    end
    meanDiffPerm(numP,:)=nanmean(MatPerm,1);
end

%observed mean difference
tmp1=nanmean(dataMat1-dataMat2,3);
meanDiffData=tmp1(Idx);  %transform matrix to linear vector

%calculate p-values
pval=nan(1,length(meanDiffData));
for comb=1:numComp
    pval(comb) =  mean(abs(meanDiffPerm(:,comb))>=abs(meanDiffData(comb)));
    if pval(comb)==0
        pval(comb)=1/numPerms; %so that p-value isn't exactly zero
    end
end
hval=pval<0.05;

AllPar.IndivCell.permUC.PVal=zeros(numparcells,numparcells);
AllPar.IndivCell.permUC.hVal=zeros(numparcells,numparcells);
AllPar.IndivCell.permUC.PVal(Idx)=pval;
AllPar.IndivCell.permUC.hVal(Idx)=hval;

[hval,~,~,pval]=fdr_bh(pval,0.01,'dep'); %benjamini hotchberg correction for multiple comparisons at q =0.01 fdr

AllPar.IndivCell.permBHC.PVal=zeros(numparcells,numparcells);
AllPar.IndivCell.permBHC.hVal=zeros(numparcells,numparcells);
AllPar.IndivCell.permBHC.PVal(Idx)=pval;
AllPar.IndivCell.permBHC.hVal(Idx)=hval;
end 
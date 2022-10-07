function [CorrMatrix]=calcMaxCrossCorr(ColorA,ColorB,timeShift)
%calculates cross-correlation between imaging data in ColorA and imaging
%data in ColorB, and time lags specified by timeShift (in frames)
numParcells=size(ColorA,1); 
CorrMatrix=nan(numParcells,numParcells); 
for par1 =1:numParcells
    for par2=1:numParcells
    %get the maximum of cross-correlation with time lag 
    if any(isnan(ColorA(par1,:))) || any(isnan(ColorB(par2,:))),continue, end 
    [tmp,~]=crosscorr(ColorA(par1,:),ColorB(par2,:),timeShift); 
    [~,I]=max(abs(tmp));
    CorrMatrix(par1,par2)=tmp(I);
    end 
end 
end 
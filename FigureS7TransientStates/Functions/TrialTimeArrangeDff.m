function[trial_dff,timeStamps]= TrialTimeArrangeDff(dff,frameTS,fsimaging,preEventWin,eventTS,postEventWin)
%This function arranges scanimage data across trials and across time based on scanimage triggers
%and sampling rate and  outputs dff and corresponding timestamps. If
%frames are dropped, it ouptuts NaNs for that timepoint 
frameNumber=1:length(frameTS);
frameTS=(frameTS(:))';
combFrames=[frameNumber', frameTS']; 
finMatrix1=nan(length(eventTS),ceil(fsimaging*(preEventWin+postEventWin)));
finMatrix2=nan(length(eventTS),ceil(fsimaging*(preEventWin+postEventWin)));
for i=1:length(eventTS)
    currEvt=eventTS(i);  
    evtarray=round(((currEvt-preEventWin):1/fsimaging:(currEvt+postEventWin)),4);%round to four decimal point
    for j=1:(length(evtarray)-1)
        tmp1=combFrames(frameTS>=evtarray(j)& frameTS<evtarray(j+1),1);
        tmp2=combFrames(frameTS>=evtarray(j)& frameTS<evtarray(j+1),2);
        if isempty(tmp1), continue; end
        if length(tmp1)==1    
            finMatrix1(i,j)=tmp1;
            finMatrix2(i,j)=tmp2;
        elseif length(tmp1)==2 && j>1 && isnan(finMatrix1(i,j-1))
            finMatrix1(i,j-1)=tmp1(1);
            finMatrix2(i,j-1)=tmp2(1);
            finMatrix1(i,j)=tmp1(2);
            finMatrix2(i,j)=tmp2(2);
        else
            finMatrix1(i,j)=tmp1(1);
            finMatrix2(i,j)=tmp2(1); 
        end
    end
    clear currEvt evtarray tmp1 tmp2        
end
%initialize the final matrix 
numUnits=size(dff,1); 
finMatrix3=nan(length(eventTS),ceil(fsimaging*(preEventWin+postEventWin)),numUnits);
%organize the flourescence values trial by trial 
for k=1:numUnits
    currDff=dff(k,:);
    for r=1:size(finMatrix1,1)
        for s=1:size(finMatrix1,2)
           tmp1=finMatrix1(r,s);
           if ~isnan(tmp1)
           finMatrix3(r,s,k)=currDff(tmp1);
           end 
        end 
    end 
    clear currDff tmp1 
end 
trial_dff=permute(finMatrix3,[2 1 3]);  %switch rows and columns 
timeStamps=permute(finMatrix2,[2 1 3]);%switch rows and columns 
timeStamps=timeStamps-eventTS; 
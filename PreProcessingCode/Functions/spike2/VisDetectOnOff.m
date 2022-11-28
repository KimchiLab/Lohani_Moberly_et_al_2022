function[StimOnTS, StimOffTS]= VisDetectOnOff(data,threshold1,threshold2,stimDur,stimITI,spike2SR)
visTrace = medfilt1(data,100); %
vis=tsmovavg(visTrace,'s',50,1); %
vis=vis-nanmean(vis);
vis=(vis/max(vis)>threshold1); %use 0.6 threshold to get ride of start and end artifacts
[vison,visoff]=squaredetect(vis,threshold2);
StimOnTS=vison*1/spike2SR;
StimOffTS=visoff*1/spike2SR;
tmp2=StimOffTS-StimOnTS; % remove artifact in stim on and stim off times if they are greater than stim duration + .2 s
tmp3=find(tmp2>(stimDur+0.1) | (tmp2 <stimDur-0.1));
if ~isempty(tmp3)
    StimOnTS(tmp3)=[];
    StimOffTS(tmp3)=[];
end
tmp4=diff(StimOnTS); 
tmp5=find(tmp4>(stimDur+stimITI+0.5)| tmp4<(stimDur+stimITI-0.5)); %remove artifact in stim on /off times if the iti time is greater or less than the said iti
if ~isempty(tmp5)
    StimOnTS(tmp5)=[];
    StimOffTS(tmp5)=[];
end
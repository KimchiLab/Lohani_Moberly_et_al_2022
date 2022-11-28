
function [h1,sCFG] =wheel_changepoints(sCFG,dataWheel)
%written by Quentin Perrenoud, Cardin Lab  
sampleRate=dataWheel.fsample;
db1WheelRot=cell2mat(dataWheel.trial(1));
db1Range = [min(db1WheelRot) max(db1WheelRot)];
db1WheelVRange = [min([db1Range(1), sCFG.sPARAM.db1WheelVRange(1)]) ...
    max([db1Range(2), sCFG.sPARAM.db1WheelVRange(2)])];

%Express wheel rotation in radian and corrects for the smoothing of phase
%transition comming from experimental smoothing
db1WheelRotRad = (((db1WheelRot - db1WheelVRange(1))/range(db1WheelVRange))*2*pi) - pi;

%Finds the area whre the wheel rotation resets
bl1Reversal = abs(diff(db1WheelRotRad)) > 1;

in1BegPt = find(diff(bl1Reversal) == 1); in1EndPt = find(diff(bl1Reversal) == -1);
%Makes sure than beginning are always followed by an end
if ~isempty(in1BegPt)
    if in1BegPt(1) > in1EndPt(1); in1EndPt(1) = []; end
    if in1BegPt(end) > in1EndPt(end); in1BegPt(end) = []; end
end

for i = 1:length(in1BegPt)
    inAntBegPt = max(in1BegPt(i) - round(0.025 * sampleRate), 1);
    [dummy, inMaxIdx] = max(db1WheelRotRad(inAntBegPt:in1BegPt(i))); 
    in1BegPt(i) = inMaxIdx - 1 + inAntBegPt; %point right after the max
    inPostEndPt = min(in1EndPt(i) + round(0.025 * sampleRate), length(db1WheelRotRad));
    [dummy, inMinIdx] = min(db1WheelRotRad(in1EndPt(i):inPostEndPt)); 
    in1EndPt(i) = inMinIdx - 1 + in1EndPt(i); %point right befor the max
end

in1MidPt = round(mean([in1BegPt; in1EndPt]));

for i = 1:length(in1BegPt)
    db1WheelRotRad(in1BegPt(i):in1MidPt(i)) = pi*sign(db1WheelRotRad(in1BegPt(i)));
    db1WheelRotRad(in1MidPt(i)+1:in1EndPt(i)) = pi*sign(db1WheelRotRad(in1EndPt(i)));
end

%Unwraps the phase and normalize by diameter to express running distance
db1WheelDist = (unwrap(db1WheelRotRad));
db1RunDistM = ((db1WheelDist - db1WheelDist(1) / pi) * sCFG.sPARAM.dbWheelDiameterM)./4; %in meters
if db1RunDistM(end) < db1RunDistM(1)
    db1RunDistM = - db1RunDistM; %Puts the signal back in the right sense if it is inverted
end


%Computes the speed by smoothing the differential of the distance
db1Win = gausswin(sampleRate*0.5); %for a smoothing of half a second
db1SpeedMpS = (conv(diff(db1RunDistM)*sampleRate, db1Win, 'same')/sum(db1Win)); %speed in m per second
db1SpeedMpS = [db1SpeedMpS db1SpeedMpS(end)];

% % plots the result for the purpose of debugging (uncomment if needed)
% figure, ax(1) = subplot(3, 1, 1); plot(db1WheelRotRad), title(strrep(chDestFolder, '_', '\_'))
% ax(2) = subplot(3, 1, 2); plot(db1RunDistM), ax(3) = subplot(3, 1, 3); plot(db1SpeedMpS)
% linkaxes(ax, 'x')

% Writes the output in CFG
sCFG.sL0PPWR.db1WheelRotRad = db1WheelRotRad;
sCFG.sL0PPWR.db1RunDist = db1RunDistM;
sCFG.sL0PPWR.db1SpeedMpS = db1SpeedMpS;
sCFG.sL0PPWR.db1TStamps = cell2mat(dataWheel.time(1));
sCFG.sL0PPWR.sampleRate = sampleRate;
sCFG.sL0PPWR.chScriptName = mfilename('fullpath');
sCFG.sL0PPWR.chTimeComputed = datestr(now);

%% identify change points 
db1SpeedMpS     = sCFG.sL0PPWR.db1SpeedMpS;
db1TStamps      = sCFG.sL0PPWR.db1TStamps;
inSampleRate    = sCFG.sL0PPWR.sampleRate;

%Sets som parameters 
dbWindowLenSec = sCFG.sPARAM.dbWindowLenSec; % set the temoporal resolution of the analysis
bl1IsNumeric = double(~isnan(db1SpeedMpS)); %This will be usefull to compute the moving SD

%Computes the sd over a moving window (this part is the most
%computationaly expensive)
inWLen = round(dbWindowLenSec*inSampleRate);
if mod(inWLen, 2) ~= 0, inWLen = inWLen + 1; end %Makes sure the window length is even
in1Win  = ones(1,inWLen);
db1MovAvrgSpeed         = conv2(db1SpeedMpS, in1Win,'same');
in1DoF                  = conv2(bl1IsNumeric, in1Win,'same');
db1MovAvSpeedSquare     = conv2(db1SpeedMpS.^2, in1Win,'same');
db1SpeedMovSD           = sqrt((db1MovAvSpeedSquare - db1MovAvrgSpeed.^2./in1DoF)./(in1DoF-1));%moving standard deviation with inWLen as the window size

%renormalizes the moving average speed
db1MovAvrgSpeed = db1MovAvrgSpeed./in1DoF;

%Computes a forward zscore corresponding to the mean of the speed in the
%window after the point expressed in units of the standard deviation of the
%speed in the window before the point. Computes a backward zscore following
%a similar principle.
db1ZForward  = [nan(1, inWLen) db1MovAvrgSpeed(inWLen + (inWLen/2) + 1:end -(inWLen/2))...
    ./db1SpeedMovSD((inWLen/2):end - (inWLen/2) - inWLen - 1) nan(1, inWLen)];
db1ZBackward  = [nan(1, inWLen) db1MovAvrgSpeed(inWLen/2 :end -(inWLen/2) - inWLen - 1)...
    ./db1SpeedMovSD(inWLen + inWLen/2:end - inWLen/2 - 1) nan(1, inWLen)];
%Threshold the moving SD to get a crude estimate of the period of mobility
%and a first estimate of transition points
% fprintf('Moving SD SD: %.4f    \r', std(db1SpeedMovSD)), return
bl1ZSDSig           = db1SpeedMovSD > 0.005;  %the threshold is empirical and might change if the moving speed is smoothed differently 
bl1ZSDSig([1 end])  = 0; %Makes sure that the recording begins and ends with periods of immobility
in1OnIdx    = find(diff(bl1ZSDSig) == 1); %first estimates of 'wheel on'
in1OffIdx   = find([false diff(bl1ZSDSig) == -1]); %first estimate 'wheel off'
%Gets a better estimate using the maximum of the moving zscore in a window
%around each points whose length is defined by window length sec
for ii = 1:length(in1OnIdx)
    if in1OnIdx(ii) == 1, continue, end %skips the beginning of the recording if it is a period of mouvement
    inBegWin    = max(in1OnIdx(ii) - inWLen/2, 1);
    inEndWind   = min(in1OnIdx(ii) + inWLen/2, length(db1SpeedMpS));
    [ignore, inMaxIdx] = max(db1ZForward(inBegWin : inEndWind));
    in1OnIdx(ii) = in1OnIdx(ii) - inWLen/2 + inMaxIdx - 1;
end
for ii = 1:length(in1OffIdx)
    if in1OffIdx(ii) == length(db1SpeedMpS), continue, end %skips the end of the recording if it is a period of mouvement
    inBegWin    = max(in1OffIdx(ii) - inWLen/2, 1);
    inEndWind   = min(in1OffIdx(ii) + inWLen/2, length(db1SpeedMpS));
    [ignore, inMaxIdx] = max(db1ZBackward(inBegWin : inEndWind));
    in1OffIdx(ii) = in1OffIdx(ii) - inWLen/2 + inMaxIdx - 1;
end

%Removes periods of movement if they are too short
in1RemIdx = find(in1OffIdx - in1OnIdx <= inWLen/2);
in1OnIdx(in1RemIdx) = [];
in1OffIdx(in1RemIdx) = [];

if ~isempty(in1OnIdx) %SL added this if in1OnIdx can be less than 0 
    in1OnIdx(find(in1OnIdx<1))=1;
    in1OffIdx(find(in1OffIdx<1))=1;
end


%Checks that the average speed during moving period is not to low 
db1MovingAvSpeed   = nan(size(in1OnIdx));
for ii = 1:length(in1OffIdx)
    db1MovingAvSpeed(ii) = mean(abs(db1SpeedMpS(in1OnIdx(ii):in1OffIdx(ii))));
end
in1RemIdx = find(db1MovingAvSpeed <= 0.02);
in1OnIdx(in1RemIdx) = [];
in1OffIdx(in1RemIdx) = [];

if sCFG.sPARAM.blDoPlot
h1=figure, subplot(2,1,1)
    lowerLim=min(sCFG.sPARAM.blDoPlotDuration);
    upperLim=max(sCFG.sPARAM.blDoPlotDuration);
    upperLim=min(upperLim,length(dataWheel.trial{1})); 
    plot(lowerLim/sampleRate:1/sampleRate:(upperLim/sampleRate),db1WheelRot(lowerLim:upperLim)),db1YL = ylim;
    hold on;
    x=in1OnIdx(in1OnIdx>=lowerLim & in1OnIdx<upperLim);
    y=in1OffIdx(in1OffIdx>=lowerLim & in1OffIdx<upperLim);
    for ii = 1:length(x)
        plot(([x(ii)/sampleRate x(ii)/sampleRate]), db1YL, 'g')
    end
    for ii = 1:length(y)
        plot(([y(ii)/sampleRate y(ii)/sampleRate]), db1YL, 'r')
    end
    xlabel('time (s)'), ylabel('Wheel Rotation)')

subplot(2,1,2) 
    plot(lowerLim/sampleRate:1/sampleRate:(upperLim/sampleRate),db1SpeedMpS(lowerLim:upperLim)),db1YL = ylim;
    hold on;
    x=in1OnIdx(in1OnIdx>=lowerLim & in1OnIdx<upperLim);
    y=in1OffIdx(in1OffIdx>=lowerLim & in1OffIdx<upperLim);
    for ii = 1:length(x)
        plot(([x(ii)/sampleRate x(ii)/sampleRate]), db1YL, 'g')
    end
    for ii = 1:length(y)
        plot(([y(ii)/sampleRate y(ii)/sampleRate]), db1YL, 'r')
    end
    xlabel('time (s)'), ylabel('Speed (m/s)')

end 

  
%Writes the output 
sCFG.sL1DWCP.in1WheelOnIdx      = in1OnIdx;
sCFG.sL1DWCP.in1WheelOffIdx     = in1OffIdx;
sCFG.sL1DWCP.db1WheelOnTStamp   = db1TStamps(in1OnIdx);
sCFG.sL1DWCP.db1WheelOffTStamp  = db1TStamps(in1OffIdx);
end 
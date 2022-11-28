%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads spike2 smr file and saves it as a mat file
% input:    
%           data_time_stamp_filename - smr file
%           channels_num             - number of channels
%           analog_rate              - sampling rate
% outputs:
%           data                     - a matrix of all channels from smr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [data] = process_spike2_smr2mat(datapath, outputpath, data_time_stamp_filename, channels_num)

display(strcat('loading in smr: ',data_time_stamp_filename));
if ~exist(fullfile(datapath, strcat(data_time_stamp_filename)),'file')
        error('Could not find file');
end
fhand = CEDS64Open(fullfile(datapath, strcat(data_time_stamp_filename)));
if fhand == -1
    error('Could not open file');
end
ichannum = min(CEDS64MaxChan(fhand),channels_num);
[~, filename] = fileparts(data_time_stamp_filename);

if ~exist(fullfile(outputpath, strcat(filename,'.mat')), 'file')
    maxTimeTicks = CEDS64ChanMaxTime(fhand,1);
    data=nan(maxTimeTicks./20+1,ichannum);
    
    % get waveform data from each channel
    for ichan=1:ichannum
        %file name, channel num, max num of points, start and end time in ticks
        [fRead,fVals,fTime] = CEDS64ReadWaveF(fhand,ichan,maxTimeTicks,0,maxTimeTicks);
        if fRead > 0
        data(1:fRead,ichan)=fVals;
        readF(ichan) = fRead;
        end
    end
%     if length(unique(readF))~=1
%         error('Not all channels have the same length!');
%     end
    save(fullfile(outputpath, strcat(filename,'.mat')),'data','-v7.3');
else
    load(fullfile(outputpath, strcat(filename,'.mat')),'data');
end
CEDS64CloseAll();



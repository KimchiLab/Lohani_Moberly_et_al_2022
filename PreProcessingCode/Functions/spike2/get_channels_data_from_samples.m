function [channels_data,wheelOn,wheelOff,h1] = get_channels_data_from_samples(data, channels,sCFG,fsspike2)

names = fieldnames(channels);
for ni = 1:length(names)
    switch names{ni}
        case 'MOVING_PHOTO_DIODE'
            channels_data.movingdiode=data(:,channels.MOVING_PHOTO_DIODE)-nanmedian(data(:,channels.MOVING_PHOTO_DIODE));
            tmid = round(length(channels_data.movingdiode)/2);
            maxval=nanmean(channels_data.movingdiode(1:1000));
            channels_data.movingdiode(channels_data.movingdiode>maxval)=0;
        case 'BLUE'
            channels_data.blue=(data(:,channels.BLUE)-nanmin(data(:,channels.BLUE))>0.5);
        case 'UV'
            channels_data.uv=(data(:,channels.UV)-nanmin(data(:,channels.UV))>0.5);
        case 'GREEN'
            channels_data.green=(data(:,channels.GREEN)-nanmin(data(:,channels.GREEN))>0.5);
        case {'FRAMETICKS' 'RED_MESO'}
            channels_data.mesoframe=(data(:,channels.FRAMETICKS)-nanmin(data(:,channels.FRAMETICKS))>0.5);
        case 'PHOTO_DIODE'
            %channels_data.diode = extract_photo_diod_signal(data(:,channels.PHOTO_DIODE), rf);
            channels_data.diode=data(:, channels.PHOTO_DIODE);
        case 'WHEEL'           
            dataWheel.fsample=fsspike2;
            wheelData=data(:,channels.WHEEL);
            wheelData=fillmissing(wheelData,'nearest');
            dataWheel.trial{1}=(wheelData(:))'; 
            dataWheel.time{1}=1/dataWheel.fsample:1/dataWheel.fsample:((size(dataWheel.trial{1},2))/dataWheel.fsample);
            [h1,sCFG] = wheel_changepoints(sCFG,dataWheel);
            channels_data.wheelspeed=sCFG.sL0PPWR.db1SpeedMpS;
            wheelOn=sCFG.sL1DWCP.db1WheelOnTStamp;
            wheelOff=sCFG.sL1DWCP.db1WheelOffTStamp;
        case 'AIR_PUFF'
            % Channel6- Airpuff
            channels_data.air_puff=(data(:,channels.AIR_PUFF)-nanmean(data(:,channels.AIR_PUFF))>1); %airpuff
        case 'PUPIL_CAMERA'
            % Channel7- get pupil camera
             channels_data.pupil=data(:,channels.PUPIL_CAMERA);
        case 'EEG'
            channels_data.EEG=data(:,channels.EEG);
            
        otherwise
            error('Unindetified channel name');
    end
end




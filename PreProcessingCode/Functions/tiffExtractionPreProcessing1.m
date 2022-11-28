function [sigsMov, R,C]= tiffExtractionPreProcessing1(tiffsPath,outputPath, TmpFile,params)
%%written by Hadas Benisty and Sweyta Lohani 2020 
%Does the first step of Mesoscopic 1P data processing by loading tiff image sequence or cxd image and separating into dfferent colors
if nargin == 4
    batchSize = params.batchSize;
    batches2load=params.batches2load;
    moveLocal=params.moveLocal; 
    resizeRatio=params.resizeRatio;
    imageFormat=params.ImageFormat;  
else
    batchSize = 3000;
    batches2load=1000;
    moveLocal=1; 
    imageFormat='tif';
    
    % for blue-uv settings
    if strcmp(params.signalsExtraction.sigs,'blueuv')
    params.signalsExtraction.sigs = 'blueuv';
    params.signalsExtraction.firstCaFrame = 1;
    params.signalsExtraction.blueuvRatio = 1;
    resizeRatio = 0.5;
    end 
    
    if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
    params.signalsExtraction.sigs = 'RCaMP_AC';
    params.signalsExtraction.blueuvRatio = 1;
    params.signalsExtraction.firstCaFrame = 1;
    resizeRatio = 0.5;
    end 
end

%% load raw tiff movies 
if exist(fullfile(outputPath, 'raw_mov.mat'), 'file')
   load(fullfile(outputPath, 'raw_mov.mat'),'sigsMov','R','C');  
else
    % load tiffs
    mkdir(outputPath);
    disp(['load tiffs from ' tiffsPath ]);
    if strcmp(imageFormat,'tif')
    mov = read_batches_tiffs(batchSize, tiffsPath, batches2load, resizeRatio,TmpFile,moveLocal);
    elseif strcmp (imageFormat,'cxd')
    mov=readCxdSingleFile(tiffsPath,resizeRatio);
    else
    error('Unknown file format'); 
    end 
    [R, ~, ~] = size(mov);
    %% separating channels
    disp('Extract signals from tiffs');
    sigsMov = extract_signals_from_mov(mov, params.signalsExtraction);
    C=size(sigsMov.blue,1)./R; 
    save(fullfile(outputPath, 'raw_mov.mat'),'sigsMov','R','C', '-v7.3');
end 

end 

 
function [dFoF,dFoF_parcells, R,C]= tiffExtractionPreProcessing2(tiffsPath,outputPath,sigsMov,R,C,params)
%%written by Hadas Benisty 2019 and Sweyta Lohani 2020
%Does the second step of Meso 1P data processing by blue/green between camera registration, detrending
%calculation of dF/F, followed by Allen atlas parcellation
%outputs uncorrected but detrended dFoF
if nargin <5
    % detrending params
    params.deterend.filtLen = 1000;
    params.deterend.filtcutoff = 0.001;
    params.deterend.method = 'FIR';
    params.angle2rotate = -180;
end

%% Registration
C=size(sigsMov.blue,1)./R;
%if the data contains red-green-uv alternating frames
if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
    disp(['Registering green and blue images for' tiffsPath]);
    %register blue and green images
    tformfile = fullfile(outputPath, 'tform_bluegreen.mat' );
    if ~exist(tformfile, 'file')
        [tform] = registerGreenBlue(sigsMov.green,sigsMov.blue,R,C,outputPath);%register two color images from two cameras
        save(tformfile, 'tform','R','C');
    else
        load(tformfile, 'tform','R','C');
    end
    sigsMov.green = transform_frames(sigsMov.green, tform, R, C);
    sigsMov.green=single(sigsMov.green);
end
%% De-trend followed by df/f
names = fieldnames(sigsMov);
ind1 = strcmp(names,'skippedframes' );
ind2 = strcmp(names,'skippedchannels' );
siginds = find(~ind1 & ~ind2);
for i = 1:length(siginds)
    sigsMov.(names{siginds(i)}) = imrotate_vectored(sigsMov.(names{siginds(i)}),R, C, params.angle2rotate);%rotate the movie so the front of the brain is at the bottom
    firstFrame.(names{siginds(i)})=(sigsMov.(names{siginds(i)})(:,1000)); %extract one frame from the raw movie
    firstFrame.(names{siginds(i)})=imadjust(mat2gray( firstFrame.(names{siginds(i)}))); %this first frame is for plotting the parcel overlay in the parcellation step, looks better than df/f for identifying anatomical landmarks
    h=figure; subplot(1,2,1),imagesc(reshape(firstFrame.(names{siginds(i)}),R,C));title('Before'); %draw a figure before removal of saturated pixels
    maxPixelValue=(max(sigsMov.(names{siginds(i)}),[],2));
    saturatedPixelIdx=find(maxPixelValue>=65535); %find saturated pixels
    sigsMov.(names{siginds(i)})(saturatedPixelIdx,:)=0;  %remove any saturated pixel
    subplot(1,2,2),imagesc(reshape(sigsMov.(names{siginds(i)})(:,1000),R,C));%draw a figure after removal of saturated pixels
    runTitle=strcat('After',' number of saturated Pixels = ',num2str(numel(saturatedPixelIdx)));
    title(runTitle); %imageName=fullfile(outputPath,strcat('RemovalSaturatedPixels',names{i})); saveas(h,imageName);
    disp(['Detrending ' names{siginds(i)}]);
    [sigsMov.(names{siginds(i)}),sigsBaseline.(names{siginds(i)})]=detrend_all_pixels(sigsMov.(names{siginds(i)}), params.deterend);
    disp(['Extracting dFoF ' names{siginds(i)}]);
    dFoF.(names{siginds(i)}) = sigsMov.(names{siginds(i)})./ sigsBaseline.(names{siginds(i)});%calculated df/f with f0 as the low pass filtered signal
end
clear mov sigsMov
%% register movies to allen atlas and parcellate into brain regions
disp(['Atlas Registration for' tiffsPath]);
load('parcells_updated121519.mat'); parcells=parcells_new;parcells_template=(mat2gray(parcells.CombinedParcells));%load new parcells
tformfile = fullfile(outputPath, 'tform_blue.mat' );
if ~exist(tformfile, 'file')
    [tform,R,C,movingPoints,fixedPoints] = get_alignment_CotrolPoints_transform_prompt(firstFrame.blue,parcells_template,R,C,parcells.movingPoints,parcells.fixedPoints); % do alignment based on control points manually selected on GUI, you can pass what template you want to use
    if strcmp(params.signalsExtraction.sigs,'RCaMP_AC') %if there is an image in green, use that for registration instead
        [tform,R,C] = get_alignment_CotrolPoints_transform_prompt(firstFrame.green,parcells_template,R,C,movingPoints,fixedPoints); % do alignment based on control points manually selected on GUI, you can pass what template you want to use
    end
    save(tformfile, 'tform','R','C');
else
    load(tformfile, 'tform','R','C');
end

names = fieldnames(dFoF);
for i = 1:length(names)
    dFoF.(names{i}) = transform_frames(dFoF.(names{i}), tform, R, C);
    firstFrame_t.((names{i}))=transform_frames(firstFrame.(names{i}), tform, R, C);
    [h1]=plot_parcell_overlay(firstFrame_t.(names{i}),R,C,1,parcells.indicators);
    if ~isempty(h1)
        imageName1=fullfile(outputPath,strcat('Parcels_',names{i}));
        saveas(h1,imageName1);
    end
end

%% parcellate data
for i = 1:length(names)
    disp(['Parcellating' names{i}]);
    dFoF_parcells.(names{i}) = pixels2rois(dFoF.(names{i}), parcells);
end


end
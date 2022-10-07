clear all; 
%% hemodynamic correction method comparison for one example dual GRABs-RCaMP mouse session or GFP mouse session or mcherry session 
Session='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\GFP\SL225_egfp\SL225_Spont_Day3_Sess2';%one example gfp mouse session SL225_Egfp_day3_sess2 (airpuff session)
outputPath='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\Figures\Hemodynamics\SL225_egfp_Spont_Day3_Sess2'; 
params.angle2rotate=90;
params.signalsExtraction.sigs='blueuv';

% Session='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice\grabAM05\imaging with 575 excitation\11232019_grabAM05_airpuffs'; %one example GRABs-RCaMP Session (airpuff session)
% outputPath='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\Figures\Hemodynamics\11232019_grabAM05_airpuffs'; 
% params.angle2rotate=-180;
% params.signalsExtraction.sigs='RCaMP_AC';
% %  

% Session='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\mCherry\SL223_mcherry\12102019_SL223_mcherry_airpuffs3'; %mcherry session 
% outputPath='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\Figures\Hemodynamics\12102019_SL223_mcherry_airpuffs3'; 
% params.angle2rotate=-180;
% params.signalsExtraction.sigs='RCaMP_AC';

%%
params.deterend.filtLen = 1000;
params.deterend.filtcutoff = 0.001;
params.deterend.method = 'FIR';
addpath(genpath('Z:\GRABS_Data\Code\meso_processing-master'));
addpath(genpath('F:\Sweyta\GRABS_Data\FinalPaperCodes\HemodynamicCorrection')); 
load(fullfile(Session,'raw_mov.mat')); %load raw output 
load('parcells_updated121519.mat'); parcells=parcells_new;
if ~exist(outputPath,'dir'), mkdir(outputPath); end 
%% Dual cmaera registration if the data contains red-green-uv alternating frames
if ~exist(fullfile(outputPath, 'dFOF_Uncorrected'),'dir')
    if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
        %register blue and green images
        tformfile = fullfile(Session, 'tform_bluegreen.mat' );
        load(tformfile, 'tform','R','C');
        sigsMov.green = transform_frames(sigsMov.green, tform, R, C);
        sigsMov.green=single(sigsMov.green);
    end
    
    %% De-trend followed by df/f
    R=256; C=256;
    names = fieldnames(sigsMov);
    ind1 = strcmp(names,'skippedframes' );
    ind2 = strcmp(names,'skippedchannels' );
    
    siginds = find(~ind1 & ~ind2);
    for i = 1:length(siginds)
        sigsMov.(names{siginds(i)}) = imrotate_vectored(sigsMov.(names{siginds(i)}),R, C, params.angle2rotate);%rotate the movie so the front of the brain is at the bottom
        maxPixelValue=(max(sigsMov.(names{siginds(i)}),[],2));
        maxFrame1.(names{siginds(i)})=maxPixelValue;%extract maximum image for plotting purpose
        saturatedPixelIdx=find(maxPixelValue>=65535); %find saturated pixels 
        sigsMov.(names{siginds(i)})(saturatedPixelIdx,:)=0;  %remove any saturated pixel       
        disp(['Detrending ' names{siginds(i)}]);
        [sigsMov.(names{siginds(i)}),sigsBaseline.(names{siginds(i)})]=detrend_all_pixels(sigsMov.(names{siginds(i)}), params.deterend);
        disp(['Extracting dFoF ' names{siginds(i)}]);
        dFoF.(names{siginds(i)}) = sigsMov.(names{siginds(i)})./ sigsBaseline.(names{siginds(i)});%calculated df/f with f0 as the low pass filtered signal
    end
    clear mov sigsBaseline sigsMov; 
    %% register movies to allen atlas and parcellate into brain regions
    disp(['Atlas Registration']);
    parcells_template=(mat2gray(parcells.CombinedParcells));%load new parcells
    tformfile = fullfile(Session, 'tform_blue.mat' );
    load(tformfile, 'tform','R','C');
    names = fieldnames(dFoF);
    for i = 1:length(names)
        dFoF.(names{i}) = transform_frames(dFoF.(names{i}), tform, R, C);
        maxFrame_t.((names{i}))=transform_frames(maxFrame1.(names{i}), tform, R, C);
        [h1]=plot_parcell_overlay(maxFrame_t.(names{i}),R,C,1,parcells.indicators);
    end
    dFoF_Raw=dFoF; clear dFoF; 
    %save(fullfile(outputPath, 'dFOF_Raw'), 'dFoF_Raw','-v7.3');
     %get dFoF parcells
     names=fieldnames(dFoF_Raw); 
    for i = 1:length(names)
        disp(['Parcellating' names{i}]);
        dFoF_Raw_parcells.(names{i}) = pixels2rois(dFoF_Raw.(names{i}), parcells);
    end
   % save(fullfile(outputPath, 'dFOF_Raw_parcells'), 'dFoF_Raw_parcells','-v7.3');
else
    load(fullfile(outputPath, 'dFOF_Raw'), 'dFoF_Raw');
    load(fullfile(outputPath, 'dFOF_Raw_parcells'), 'dFoF_Raw_parcells');
end

%% calculate df/f after pixelwise regression with uv
if ~exist(fullfile(outputPath, 'dFOF_PixUV'), 'dir')
    names=fieldnames(dFoF_Raw);
    for i = 1:length(names)
        if  ~strcmp(names{i},'uv')
            len = min(size(dFoF_Raw.(names{i}),2), size(dFoF_Raw.uv,2));
            npix = size(dFoF_Raw.(names{i}),1);
            regsig = zeros(npix,len);
            for ipix = 1:npix
                PixxTime_bl=dFoF_Raw.(names{i})(ipix,1:len);
                PixxTime_uv=dFoF_Raw.uv(ipix,1:len);
                if any (PixxTime_bl>0)
                    al_coeff = regress(PixxTime_bl',PixxTime_uv');
                    regsig(ipix,:) = PixxTime_bl-al_coeff*PixxTime_uv;
                end
            end
            dFoF_PixUV.(names{i})=regsig;
        end
    end
    %save(fullfile(outputPath, 'dFOF_PixUV'), 'dFoF_PixUV','-v7.3');
    clear regsig;
    %get dFoF parcells
    names=fieldnames(dFoF_PixUV); 
    for i = 1:length(names)
        disp(['Parcellating' names{i}]);
        dFoF_PixUV_parcells.(names{i}) = pixels2rois(dFoF_PixUV.(names{i}), parcells);
    end
    %save(fullfile(outputPath, 'dFOF_PixUV_parcells'), 'dFoF_PixUV_parcells','-v7.3');
else
    load(fullfile(outputPath, 'dFOF_PixUV'), 'dFoF_PixUV');
    load(fullfile(outputPath, 'dFOF_PixUV_parcells'), 'dFoF_PixUV_parcells');
end


%calculate df/f after spatial regression with svd with uv (just load these because the outputs already exist) 
%load(fullfile(Session,'final_dFoF.mat')); %load dfof output 
load(fullfile(Session,'final_dFoF_parcels.mat')); %load dfof output 

%% calculate correlation matrix for the whole session 
% make correlation matrices using parcell indices in the left hemisphere only 
Idx.visual=2:2:16; %visual parcells
Idx.PosteriorParietal=32;% PPC
Idx.RSL_AgL=[18,20];%retrosplenial cortex
Idx.Somat=34:2:48;%somatosensory areas
Idx.FrontalMotor=50:2:52;%frontal-moto area
Idx.auditory=[28,30];%auditory areas
CombinedParcellIdx=[Idx.visual,Idx.RSL_AgL,Idx.auditory,Idx.PosteriorParietal,Idx.Somat,Idx.FrontalMotor];
V1Idx=1; S1bIdx=14; M2Idx=23; %updated parcell idx for left V1, S1, and M2
parcellnames=parcells.names(CombinedParcellIdx);

if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
    names={'blue','green'}; 
else
    names={'blue'}; 
end 


for i=1:length(names)
%uncorrected correlation matrix 
corrMatrix.raw.(names{i})=corr(dFoF_Raw_parcells.(names{i})(CombinedParcellIdx,:)');

%pixelwise uv regression correlation matrix
corrMatrix.pixUV.(names{i})=corr(dFoF_PixUV_parcells.(names{i})(CombinedParcellIdx,:)');

%spatial SVD correlation matrix
corrMatrix.svdUV.(names{i})=corr(dFoF_parcells.(names{i})(CombinedParcellIdx,:)');

%make correlation matrix figures 
numparcells=length(CombinedParcellIdx);  tmp=ones(numparcells,numparcells);
tmp=tril(tmp-diag(diag(tmp)));
[Indx]=find(tmp==0);

corrMatrix.raw.(names{i})(Indx)=0;
corrMatrix.pixUV.(names{i})(Indx)=0;
corrMatrix.svdUV.(names{i})(Indx)=0;
cmap = [1 1 1;parula(200)];
h1= figure;
subplot(1,3,1); imagesc(corrMatrix.raw.(names{i})); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:)); 
caxis([0 1]); title(strcat((names{i}),'Uncorrected-DF/F')); box off;
xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');

subplot(1,3,2); imagesc(corrMatrix.pixUV.(names{i})); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:)); 
caxis([0 1]); title(strcat((names{i}),'PixelwiseUVRegression-DF/F')); box off;
xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');

subplot(1,3,3); imagesc(corrMatrix.svdUV.(names{i})); xticks(1:length(parcellnames));yticks(1:length(parcellnames)); xticklabels(parcellnames); yticklabels(parcellnames(:)); 
caxis([0 1]); title(strcat((names{i}),'SpatialSVDUV-DF/F')); box off;
xtickangle(90); colormap(cmap); colorbar;set(gca,'TickDir','out');
saveas(h1,fullfile(outputPath,strcat('CorrMatrixSession-',names{i}))); 
end 
%% plot remaining variance after correction and some raw traces 
% plot traces for each color 
load(fullfile(Session, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps');
timeframes=1500:3000;
imaging_timestamps=timestamps.timaging(timeframes);
airpuff_times=timestamps.airpuff(timestamps.airpuff>imaging_timestamps(1) & timestamps.airpuff<imaging_timestamps(end));
for i=1:length(names)
h2=figure;ax(1)=subplot(3,1,1);plot(imaging_timestamps,dFoF_Raw_parcells.(names{i})(2,timeframes));  title(strcat('Uncorrected-DF/F-V1-',(names{i})));  box off;axis off;ylim([-0.04 0.02]); 
hold on;ylimits=ylim; 
for tt=1:length(airpuff_times)
    plot([airpuff_times(tt),airpuff_times(tt)], ylimits, 'k');
end
ax(2)=subplot(3,1,2);plot(imaging_timestamps,dFoF_PixUV_parcells.(names{i})(2,timeframes)); title('PixelwiseRegress-V1'); box off; axis off; ylim([-0.04 0.02]); 
hold on;ylimits=ylim; 
for tt=1:length(airpuff_times)
    plot([airpuff_times(tt),airpuff_times(tt)], ylimits, 'k');
end
ax(3)=subplot(3,1,3);plot(imaging_timestamps,dFoF_parcells.(names{i})(2,timeframes));  title('SpatialSVDUV-V1');box off; ylim([-0.04 0.02]); xlabel('Frames'); 
hold on;ylimits=ylim; 
for tt=1:length(airpuff_times)
    plot([airpuff_times(tt),airpuff_times(tt)], ylimits, 'k');
end
linkaxes(ax,'xy'); 
saveas(h2,fullfile(outputPath,strcat('RawTraces-',names{i}))); 

% % get a whole brain variance plot parcel-wise and pixel-wise-- remaining variance as fraction of raw variance 
% varianceParcelSVDUV=var(dFoF_parcells.(names{i}),0,2); 
% varianceParcelPixUV=var(dFoF_PixUV_parcells.(names{i}),0,2); 
% rawvarianceParcel=var(dFoF_Raw_parcells.(names{i}),0,2); 
% variancePixSVDUV=var(dFoF.(names{i}),0,2); 
% variancePixPixUV=var(dFoF_PixUV.(names{i}),0,2); 
% rawvariancePix=var(dFoF_Raw.(names{i}),0,2); 
% 
% fractionVarianceSVDUV_parcel=varianceParcelSVDUV./rawvarianceParcel; 
% h3=figure;subplot(1,2,1),bar(fractionVarianceSVDUV_parcel(CombinedParcellIdx));ylim([0 1]); ylabel('FractionRemainingVariance-SpatialSVDUV'); title(strcat('Parcells-',(names{i}))); 
% xticks(1:length(parcellnames)); xticklabels(parcellnames); xtickangle(90); 
% 
% fractionVarianceSVDUV_pix=variancePixSVDUV./rawvariancePix; 
% load('brainMask.mat','mask');
% maskinds = find(mask==0);
% varianceBrain=fractionVarianceSVDUV_pix; varianceBrain(maskinds)=nan; 
% varianceBrain=reshape(varianceBrain,R,C);
% subplot(1,2,2); imagesc(varianceBrain,'AlphaData',~isnan(varianceBrain));colormap('parula'); box off; axis off; caxis([0 0.8]);  title(strcat('FractionRemainingVariance-SpatialSVDUV-',(names{i}))); colorbar; 
% 
% saveas(h3,fullfile(outputPath,strcat('FractionRemainingVariance-SpatialSVDUV-',names{i}))); 
% 
% fractionVariancePixUV_parcel=varianceParcelPixUV./rawvarianceParcel; 
% h4=figure;subplot(1,2,1),bar(fractionVariancePixUV_parcel(CombinedParcellIdx));ylim([0 1]); ylabel('FractionRemainingVariance-PixelwiseRegressionUV'); title(strcat('Parcells-',(names{i}))); 
% xticks(1:length(parcellnames)); xticklabels(parcellnames); xtickangle(90); 
% 
% fractionVariancePixUV_pix=variancePixPixUV./rawvariancePix; 
% varianceBrain=fractionVariancePixUV_pix; varianceBrain(maskinds)=nan; 
% varianceBrain=reshape(varianceBrain,R,C); 
% subplot(1,2,2); imagesc(varianceBrain,'AlphaData',~isnan(varianceBrain));colormap('parula'); box off; axis off; caxis([0 0.8]); title(strcat('FractionRemainingVariance-PixelwiseRegressionUV-',(names{i}))); colorbar; 
% saveas(h4,fullfile(outputPath,strcat('FractionRemainingVariance-PixelwiseRegressionUV-',names{i}))); 
end 


%% plot averaged response to airpuff trials 

eventTS=timestamps.airpuff(:)';preEventWin=2; postEventWin=5;baselineWin=[-2 0]; fsimaging=10; timaging=timestamps.timaging; V1Idx=2;
baselineFrames=(((baselineWin(1)+preEventWin)*fsimaging)+1):((baselineWin(2)+preEventWin)*fsimaging);

for i=1:length(names)
    h5=figure; title(strcat('Airpuff-',names{i}));hold on; 
    %uncorrected raw data
    [Data_NonNorm_V1]= TrialTimeArrangeDff(dFoF_Raw_parcells.(names{i})(V1Idx,:),timaging,fsimaging,preEventWin,eventTS,postEventWin);
    Data_NonNorm_V1=fillmissing(Data_NonNorm_V1,'linear',1,'EndValues','nearest');% fill nan values
    baseMean=nanmean(Data_NonNorm_V1(baselineFrames,:)); baseStd=nanstd(Data_NonNorm_V1(baselineFrames,:));
    Data_ZNorm_V1= (Data_NonNorm_V1-baseMean)./baseStd;
    ax1= plot((-preEventWin+1/fsimaging):1/fsimaging:(postEventWin),nanmean(Data_ZNorm_V1,2)); 
    
    %pixelwise uv regression 
    [Data_NonNorm_V1]= TrialTimeArrangeDff(dFoF_PixUV_parcells.(names{i})(V1Idx,:),timaging,fsimaging,preEventWin,eventTS,postEventWin);
    Data_NonNorm_V1=fillmissing(Data_NonNorm_V1,'linear',1,'EndValues','nearest');% fill nan values
    baseMean=nanmean(Data_NonNorm_V1(baselineFrames,:)); baseStd=nanstd(Data_NonNorm_V1(baselineFrames,:));
    Data_ZNorm_V1= (Data_NonNorm_V1-baseMean)./baseStd;
    ax2=plot((-preEventWin+1/fsimaging):1/fsimaging:(postEventWin),nanmean(Data_ZNorm_V1,2)); 
    
    %sptial svd regression with uv 
    [Data_NonNorm_V1]= TrialTimeArrangeDff(dFoF_parcells.(names{i})(V1Idx,:),timaging,fsimaging,preEventWin,eventTS,postEventWin);
    Data_NonNorm_V1=fillmissing(Data_NonNorm_V1,'linear',1,'EndValues','nearest');% fill nan values
    baseMean=nanmean(Data_NonNorm_V1(baselineFrames,:)); baseStd=nanstd(Data_NonNorm_V1(baselineFrames,:));
    Data_ZNorm_V1= (Data_NonNorm_V1-baseMean)./baseStd;
    ax3= plot((-preEventWin+1/fsimaging):1/fsimaging:(postEventWin),nanmean(Data_ZNorm_V1,2)); 
    
    legend([ax1,ax2,ax3],{'Uncorrected','Pixelwise-UV','SpatialSVD-UV'}); ylabel('ZScore'); xlabel('Time (s)'); hold off; 
    saveas(h5,fullfile(outputPath,strcat('AirpuffResponse-',names{i}))); 
end
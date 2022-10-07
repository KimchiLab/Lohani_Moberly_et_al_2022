%% load file
load('F:\Sweyta\GRABS_Data\BackScatter\grab09_backscatter.mat');
outputPath='Z:\GRABS_Data\Analyzed_SVDMethodPatch14\BackScatterOutput';
addpath(genpath('Z:\GRABS_Data\Code\meso_processing-master'));
addpath(genpath('F:\Sweyta\GRABS_Data\BackScatter\')); 
%% parameters
params.angle2rotate=-180;
params.deterend.filtLen = 1000;
params.deterend.filtcutoff = 0.001;
params.deterend.method = 'FIR';

%% De-trend followed by df/f
R=256; C=256;
names = fieldnames(sigsMov);
ind1 = strcmp(names,'skippedframes' );
ind2 = strcmp(names,'skippedchannels' );

siginds = find(~ind1 & ~ind2);

for i = 1:length(siginds)
    sigsMov.(names{siginds(i)})=sigsMov.(names{siginds(i)})(:,1:9949); %remove zeros appended at the end for this specific file
    sigsMov.(names{siginds(i)}) = imrotate_vectored(sigsMov.(names{siginds(i)}),R, C, params.angle2rotate);%rotate the movie so the front of the brain is at the bottom
    maxPixelValue=(max(sigsMov.(names{siginds(i)}),[],2));
    maxFrame1.(names{siginds(i)})=maxPixelValue;%extract maximum image for plotting purpose
    disp(['Detrending ' names{siginds(i)}]);
    [sigsMov.(names{siginds(i)}),sigsBaseline.(names{siginds(i)})]=detrend_all_pixels(sigsMov.(names{siginds(i)}), params.deterend);
    disp(['Extracting dFoF ' names{siginds(i)}]);
    dFoF.(names{siginds(i)}) = sigsMov.(names{siginds(i)})./ sigsBaseline.(names{siginds(i)});%calculated df/f with f0 as the low pass filtered signal
end
clear mov sigsMov sigsBaseline

%% register movies to allen atlas and parcellate into brain regions
disp(['Atlas Registration']);
load('parcells_updated121519.mat'); parcells=parcells_new;parcells_template=(mat2gray(parcells.CombinedParcells));%load new parcells
tformfile = fullfile(outputPath, 'tform_blue.mat' );
if ~exist(tformfile, 'file')
    [tform,R,C] = get_alignment_CotrolPoints_transform_prompt(maxFrame1.blue,parcells_template,R,C,parcells.movingPoints,parcells.fixedPoints); % do alignment based on control points manually selected on GUI, you can pass what template you want to use
    save(tformfile, 'tform','R','C');
else
    load(tformfile, 'tform','R','C');
end

names = fieldnames(dFoF);
for i = 1:length(names)
    dFoF.(names{i}) = transform_frames(dFoF.(names{i}), tform, R, C);
    maxFrame_t.((names{i}))=transform_frames(maxFrame1.(names{i}), tform, R, C);
    [h1]=plot_parcell_overlay(maxFrame_t.(names{i}),R,C,1,parcells.indicators);
    if ~isempty(h1)
        imageName1=fullfile(outputPath,strcat('Parcels_',names{i}));
        saveas(h1,imageName1);
    end
end
save(fullfile(outputPath, 'dFOF_Uncorrected'), 'dFoF','-v7.3');

%% first run all approaches with gaussian smoothing of hemodynamics channels such as uv, BS1, and BS2

%smooth uv and reflectance channesl 
dFoF_sm.uv=smoothdata(dFoF.uv,2,'gaussian',5); % smooth/filter to remove high frequency noise from uv 
dFoF_sm.BS1=smoothdata(dFoF.BS1,2,'gaussian',5); % smooth/filter to remove high frequency noise from backscatter channel 1 
dFoF_sm.BS2=smoothdata(dFoF.BS2,2,'gaussian',5); % smooth/filter to remove high frequency noise from backscatter channel 1 

%% Method 1: Hadas-Boris spatial SVD method
% load brain mask
load('brainMask.mat','mask');
maskinds = find(mask);
naninds=find(isnan(dFoF.blue(:,1)));
for i=1:length(naninds)
    [naninds_sub(i,1),naninds_sub(i,2)]=ind2sub([R,C],naninds(i));
    mask(naninds_sub(i,1),naninds_sub(i,2))=0;
end
%ensure that all fields have the same number of elements
names = fieldnames(dFoF);
minsize=size(dFoF.uv,2);
for i =1:length(names)
    minsize=min(size(dFoF.(names{i}),2),minsize);
end
%% Method 1: do SVD against uv flourescence channel
%smooth uv first
dFoF_SVD_UV.uv=dFoF_sm.uv;
dFoF_SVD_UV.uv=dFoF_SVD_UV.uv(:,1:minsize);
%run spatial svd with a patch size of 14 and smoothed uv channel
for i = 1:length(names)-2 %exclude reflectance cahnnels
    if  ~strcmp(names{i},'uv')
        disp(strcat('SVD correction processing',names{i}));
        dFoF_SVD_UV.(names{i})=dFoF.(names{i})(:,1:minsize);
        [dFoF_SVD_UV.(names{i}),~ , ~] = spatial_regression(dFoF_SVD_UV.(names{i}), dFoF_SVD_UV.uv, mask, 14, 0);
    end
end

save(fullfile(outputPath, 'dFOF_SVD_UV'), 'dFoF_SVD_UV','-v7.3');
clear dFoF_SVD_UV; 
%% Method 2: do SVD against green reflectance channel
dFoF_SVD_BS1.BS1=dFoF_sm.BS1;
dFoF_SVD_BS1.BS1=dFoF_SVD_BS1.BS1(:,1:minsize);
%run spatial svd with a patch size of 14 and smoothed uv channel
for i = 1:length(names)-2 %exclude reflectance cahnnels
    if  ~strcmp(names{i},'BS1')
        disp(strcat('SVD correction processing using BS1',names{i}));
        dFoF_SVD_BS1.(names{i})=dFoF.(names{i})(:,1:minsize);
        [dFoF_SVD_BS1.(names{i}),~ , ~] = spatial_regression(dFoF_SVD_BS1.(names{i}), dFoF_SVD_BS1.BS1, mask, 14, 0);
    end
end

save(fullfile(outputPath, 'dFOF_SVD_BS1'), 'dFoF_SVD_BS1','-v7.3');
clear dFoF_SVD_BS1; 
%% Method 3: pixelwise regression against uv
dFoF_PixUV.uv=dFoF_sm.uv; 
 for i = 1:length(names)-2 %exclude backscatter channels 
     if  ~strcmp(names{i},'uv')
          len = min(size(dFoF.(names{i}),2), size(dFoF_PixUV.uv,2));
          npix = size(dFoF.(names{i}),1);          
          regsig = zeros(npix,len);
          for ipix = 1:npix
              PixxTime_bl=dFoF.(names{i})(ipix,1:len);
              PixxTime_uv=dFoF_PixUV.uv(ipix,1:len);
              if any (PixxTime_bl>0)
              al_coeff = regress(PixxTime_bl',PixxTime_uv');
              regsig(ipix,:) = PixxTime_bl-al_coeff*PixxTime_uv; 
              end 
          end
          dFoF_PixUV.(names{i})=regsig; 
     end 
 end 
 save(fullfile(outputPath, 'dFOF_PixUV'), 'dFoF_PixUV','-v7.3');
clear dFoF_PixUV regsig; 

%% Method 4: pixelwise regression against 530 backscatter
dFoF_PixBS1.BS1=dFoF_sm.BS1; 
  for i = 1:length(names)-2 %exclude backscatter channels 
          len = min(size(dFoF.(names{i}),2), size(dFoF_PixBS1.BS1,2));
          npix = size(dFoF.(names{i}),1);          
          regsig = zeros(npix,len);
           S1.(names{i})=nan(1,npix); 
          for ipix = 1:npix
              PixxTime_bl=dFoF.(names{i})(ipix,1:len);
              PixxTime_bs1=dFoF_PixBS1.BS1(ipix,1:len);
              if any (PixxTime_bl>0)&& any(~isnan(PixxTime_bs1))
              coeff = regress(PixxTime_bl',PixxTime_bs1');
              regsig(ipix,:) = PixxTime_bl-coeff*PixxTime_bs1; 
               %save S1 acoefficients to later compare with the Beer Lambert model 
              S1.(names{i})(ipix)=coeff;
              end          
          end
          dFoF_PixBS1.(names{i})=regsig; 
 end 
 save(fullfile(outputPath, 'dFOF_PixBS1'), 'dFoF_PixBS1','S1','-v7.3');
clear dFoF_PixBS1 regsig; 

%% Method 5: pixelwise regression against 530 and 620 backscatter
dFoF_PixBS12.BS1=dFoF_sm.BS1; 
dFoF_PixBS12.BS2=dFoF_sm.BS2; 
 for i = 1:length(names)-2 %exclude backscatter channels 
          len = min([size(dFoF.(names{i}),2), size(dFoF_PixBS12.BS1,2),size(dFoF_PixBS12.BS2,2)]);
          npix = size(dFoF.(names{i}),1);          
          regsig = zeros(npix,len);
          S1.(names{i})=nan(1,npix); 
          S2.(names{i})=nan(1,npix); 
          for ipix = 1:npix
              PixxTime_bl=dFoF.(names{i})(ipix,1:len);
              PixxTime_bs1=dFoF_PixBS12.BS1(ipix,1:len);
              PixxTime_bs2=dFoF_PixBS12.BS2(ipix,1:len);
              if any (PixxTime_bl>0) && any(~isnan(PixxTime_bs1)) && any (~isnan(PixxTime_bs2))
              coeff = regress(PixxTime_bl',[PixxTime_bs1',PixxTime_bs2']);
              regsig(ipix,:) = PixxTime_bl-coeff(1)*PixxTime_bs1-coeff(2)*PixxTime_bs2; 
              %save S1 and S2 coefficients to later compare with the Beer Lambert model 
              S1.(names{i})(ipix)=coeff(1); 
              S2.(names{i})(ipix)=coeff(2);
              end 
          end
          dFoF_PixBS12.(names{i})=regsig; 
 end 
 save(fullfile(outputPath, 'dFOF_PixBS12'), 'dFoF_PixBS12','S1','S2','-v7.3');
clear dFoF_PixBS12 regsig; 


%% Method 6: simplified beer lambert model (Waters 2020, similar to Hillman 2016)
%wavelengths used
lambda_EX.blue=470; lambda_EX.uv=395; lambda_EX.green=575; %excitation wavelengths 
lambda_EM.blue=520; lambda_EM.uv=520; lambda_EM.green=593; %emission wavelengths 
lambda_R1=530;%backscatter 1 wavelength 
lambda_R2=625;%backscatter 2 wavelength 

%extract the dF/F corrected based on the linear equation dfc/Fc=dF/F-S1*dI1/I1-S2*dI2/I2;
names = fieldnames(dFoF);
dFoF_BLam.BS1=dFoF_sm.BS1; 
dFoF_BLam.BS2=dFoF_sm.BS2; 
for i = 1:length(names)-2 %exclude reflectance cahnnels
        [S1.(names{i}),S2.(names{i})]=computeBeerLamberCoeffs(lambda_EX.(names{i}),lambda_EM.(names{i}),lambda_R1,lambda_R2);%computer S1 and S2 coefficients 
        dFoF_BLam.(names{i}) = dFoF.(names{i}) - S1.(names{i})*dFoF_BLam.BS1 - S2.(names{i})*dFoF_BLam.BS2; 
end

save(fullfile(outputPath, 'dFOF_BLam'), 'dFoF_BLam','S1','S2','-v7.3');

%% repeat all above approaches but without smoothing of uv and reflectance channels 
%% Method 1: Hadas-Boris spatial SVD method
% load brain mask
load('brainMask.mat','mask');
maskinds = find(mask);
naninds=find(isnan(dFoF.blue(:,1)));
for i=1:length(naninds)
    [naninds_sub(i,1),naninds_sub(i,2)]=ind2sub([R,C],naninds(i));
    mask(naninds_sub(i,1),naninds_sub(i,2))=0;
end
%ensure that all fields have the same number of elements
names = fieldnames(dFoF);
minsize=size(dFoF.uv,2);
for i =1:length(names)
    minsize=min(size(dFoF.(names{i}),2),minsize);
end
%% Method1: do SVD against uv flourescence channel
%smooth uv first
dFoF_SVD_UV.uv=dFoF.uv;
dFoF_SVD_UV.uv=dFoF_SVD_UV.uv(:,1:minsize);
%run spatial svd with a patch size of 14 and smoothed uv channel
for i = 1:length(names)-2 %exclude reflectance cahnnels
    if  ~strcmp(names{i},'uv')
        disp(strcat('SVD correction processing',names{i}));
        dFoF_SVD_UV.(names{i})=dFoF.(names{i})(:,1:minsize);
        [dFoF_SVD_UV.(names{i}),~ , ~] = spatial_regression(dFoF_SVD_UV.(names{i}), dFoF_SVD_UV.uv, mask, 14, 0);
    end
end

save(fullfile(outputPath, 'dFOF_SVD_UV_nosmooth'), 'dFoF_SVD_UV','-v7.3');
clear dFoF_SVD_UV; 
%% Method 2: do SVD against green reflectance channel
dFoF_SVD_BS1.BS1=dFoF.BS1;
dFoF_SVD_BS1.BS1=dFoF_SVD_BS1.BS1(:,1:minsize);
%run spatial svd with a patch size of 14 and smoothed uv channel
for i = 1:length(names)-2 %exclude reflectance cahnnels
    if  ~strcmp(names{i},'BS1')
        disp(strcat('SVD correction processing using BS1',names{i}));
        dFoF_SVD_BS1.(names{i})=dFoF.(names{i})(:,1:minsize);
        [dFoF_SVD_BS1.(names{i}),~ , ~] = spatial_regression(dFoF_SVD_BS1.(names{i}), dFoF_SVD_BS1.BS1, mask, 14, 0);
    end
end

save(fullfile(outputPath, 'dFOF_SVD_BS1_nosmooth'), 'dFoF_SVD_BS1','-v7.3');
clear dFoF_SVD_BS1; 
%% Method 3: pixelwise regression against uv
dFoF_PixUV.uv=dFoF.uv; 
 for i = 1:length(names)-2 %exclude backscatter channels 
     if  ~strcmp(names{i},'uv')
          len = min(size(dFoF.(names{i}),2), size(dFoF_PixUV.uv,2));
          npix = size(dFoF.(names{i}),1);          
          regsig = zeros(npix,len);
          for ipix = 1:npix
              PixxTime_bl=dFoF.(names{i})(ipix,1:len);
              PixxTime_uv=dFoF_PixUV.uv(ipix,1:len);
              if any (PixxTime_bl>0)
              al_coeff = regress(PixxTime_bl',PixxTime_uv');
              regsig(ipix,:) = PixxTime_bl-al_coeff*PixxTime_uv; 
              end 
          end
          dFoF_PixUV.(names{i})=regsig; 
     end 
 end 
 save(fullfile(outputPath, 'dFOF_PixUV_nosmooth'), 'dFoF_PixUV','-v7.3');
clear dFoF_PixUV regsig; 

%% Method 4: pixelwise regression against 530 backscatter
dFoF_PixBS1.BS1=dFoF.BS1; 
 for i = 1:length(names)-2 %exclude backscatter channels 
          len = min(size(dFoF.(names{i}),2), size(dFoF_PixBS1.BS1,2));
          npix = size(dFoF.(names{i}),1);          
          regsig = zeros(npix,len);
           S1.(names{i})=nan(1,npix); 
          for ipix = 1:npix
              PixxTime_bl=dFoF.(names{i})(ipix,1:len);
              PixxTime_bs1=dFoF_PixBS1.BS1(ipix,1:len);
              if any (PixxTime_bl>0)&& any(~isnan(PixxTime_bs1))
              coeff = regress(PixxTime_bl',PixxTime_bs1');
              regsig(ipix,:) = PixxTime_bl-coeff*PixxTime_bs1; 
              %save S1 acoefficients to later compare with the Beer Lambert model 
              S1.(names{i})(ipix)=coeff; 
              end 
          end
          dFoF_PixBS1.(names{i})=regsig; 
 end 
 save(fullfile(outputPath, 'dFOF_PixBS1_nosmooth'), 'dFoF_PixBS1','S1','-v7.3');
clear dFoF_PixBS1 regsig; 

%% Method 5: pixelwise regression against 530 and 620 backscatter
dFoF_PixBS12.BS1=dFoF.BS1; 
dFoF_PixBS12.BS2=dFoF.BS2; 
 for i = 1:length(names)-2 %exclude backscatter channels 
          len = min([size(dFoF.(names{i}),2), size(dFoF_PixBS12.BS1,2),size(dFoF_PixBS12.BS2,2)]);
          npix = size(dFoF.(names{i}),1);          
          regsig = zeros(npix,len);
          S1.(names{i})=nan(1,npix); 
          S2.(names{i})=nan(1,npix); 
          for ipix = 1:npix
              PixxTime_bl=dFoF.(names{i})(ipix,1:len);
              PixxTime_bs1=dFoF_PixBS12.BS1(ipix,1:len);
              PixxTime_bs2=dFoF_PixBS12.BS2(ipix,1:len);
              if any (PixxTime_bl>0) && any(~isnan(PixxTime_bs1)) && any (~isnan(PixxTime_bs2))
              coeff = regress(PixxTime_bl',[PixxTime_bs1',PixxTime_bs2']);
              regsig(ipix,:) = PixxTime_bl-coeff(1)*PixxTime_bs1-coeff(2)*PixxTime_bs2; 
              %save S1 and S2 coefficients to later compare with the Beer Lambert model 
              S1.(names{i})(ipix)=coeff(1); 
              S2.(names{i})(ipix)=coeff(2);
              end 
          end
          dFoF_PixBS12.(names{i})=regsig; 
 end 
 save(fullfile(outputPath, 'dFOF_PixBS12_nosmooth'), 'dFoF_PixBS12','S1','S2','-v7.3');
clear dFoF_PixBS12 regsig; 


%% Method 6: simplified beer lambert model (Waters 2020, similar to Hillman 2016)
%wavelengths used
lambda_EX.blue=470; lambda_EX.uv=395; lambda_EX.green=575; %excitation wavelengths 
lambda_EM.blue=520; lambda_EM.uv=520; lambda_EM.green=593; %emission wavelengths 
lambda_R1=530;%backscatter 1 wavelength 
lambda_R2=625;%backscatter 2 wavelength 

%extract the dF/F corrected based on the linear equation dfc/Fc=dF/F-S1*dI1/I1-S2*dI2/I2;
names = fieldnames(dFoF);
dFoF_BLam.BS1=dFoF.BS1; 
dFoF_BLam.BS2=dFoF.BS2; 
for i = 1:length(names)-2 %exclude reflectance cahnnels
        [S1.(names{i}),S2.(names{i})]=computeBeerLamberCoeffs(lambda_EX.(names{i}),lambda_EM.(names{i}),lambda_R1,lambda_R2);%computer S1 and S2 coefficients 
        dFoF_BLam.(names{i}) = dFoF.(names{i}) - S1.(names{i})*dFoF_BLam.BS1 - S2.(names{i})*dFoF_BLam.BS2; 
end

save(fullfile(outputPath, 'dFOF_BLam_nosmooth'), 'dFoF_BLam','S1','S2','-v7.3');

%% Method 7: do SVD against red reflectance channel
dFoF_SVD_BS2.BS2=dFoF.BS2;
dFoF_SVD_BS2.BS2=dFoF_SVD_BS2.BS2(:,1:minsize);
%run spatial svd with a patch size of 14 and smoothed uv channel
for i = [1,3] % for blue and green channel only 
    if  ~strcmp(names{i},'BS2')
        disp(strcat('SVD correction processing using BS2',names{i}));
        dFoF_SVD_BS2.(names{i})=dFoF.(names{i})(:,1:minsize);
        [dFoF_SVD_BS2.(names{i}),~ , ~] = spatial_regression(dFoF_SVD_BS2.(names{i}), dFoF_SVD_BS2.BS2, mask, 14, 0);
    end
end

save(fullfile(outputPath, 'dFOF_SVD_BS2_nosmooth'), 'dFoF_SVD_BS2','-v7.3');
clear dFoF_SVD_BS2; 

%% Method 8: pixelwise regression against 625 backscatter
dFoF_PixBS2.BS2=dFoF.BS2; 
 for i = 1:length(names)-2 %exclude backscatter channels 
          len = min(size(dFoF.(names{i}),2), size(dFoF_PixBS2.BS2,2));
          npix = size(dFoF.(names{i}),1);          
          regsig = zeros(npix,len);
           S1.(names{i})=nan(1,npix); 
          for ipix = 1:npix
              PixxTime_bl=dFoF.(names{i})(ipix,1:len);
              PixxTime_bs1=dFoF_PixBS2.BS2(ipix,1:len);
              if any (PixxTime_bl>0)&& any(~isnan(PixxTime_bs1))
              coeff = regress(PixxTime_bl',PixxTime_bs1');
              regsig(ipix,:) = PixxTime_bl-coeff*PixxTime_bs1; 
              %save S1 acoefficients to later compare with the Beer Lambert model 
              S1.(names{i})(ipix)=coeff; 
              end 
          end
          dFoF_PixBS2.(names{i})=regsig; 
 end 
 save(fullfile(outputPath, 'dFOF_PixBS2_nosmooth'), 'dFoF_PixBS2','S1','-v7.3');
clear dFoF_PixBS2 regsig; 


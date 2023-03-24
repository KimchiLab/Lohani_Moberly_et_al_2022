function [dFoF,dFoF_parcells, R,C]= RegressionProcessing(dFoF,R,C,params,outputPath)
%%written by Hadas Benisty 2019 and edited/modified by Sweyta Lohani
%Does the third step of Meso 1P data processing by regressing out uv using
%spatial regression method

%% Spatial SVD method for each color pair
% load brain mask
load('brainMask.mat','mask');
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
dFoF.uv=dFoF.uv(:,1:minsize);
%run spatial svd with a patch size of 14
for i = 1:length(names)
    if  ~strcmp(names{i},'uv')
        if exist(fullfile(outputPath, strcat('SVDOut',names{i},'.mat')), 'file') % if the output already exists load that % if the output already exists load that
            disp(strcat('LoadingSVDOutput ',names{i}));
            load(fullfile(outputPath, strcat('SVDOut',names{i})),'output_sig');
            dFoF.(names{i})=output_sig;clear output_sig;
        else
            disp(strcat('SVD correction processing',names{i}));
            dFoF.(names{i})=dFoF.(names{i})(:,1:minsize);
            [output_sig, sig_blue, sig_uv] = spatial_regression(dFoF.(names{i}), dFoF.uv, mask, params.patchSize, 0);
            save(fullfile(outputPath, strcat('SVDOut',names{i})),'output_sig', 'sig_blue','sig_uv','-v7.3' );
            dFoF.(names{i})=output_sig;clear output_sig;
        end
    end
end
%% parcellate data
load('parcells_updated121519.mat'); parcells=parcells_new;%load parcells 
for i = 1:length(names)
    disp(['Parcellating' names{i}]);
    dFoF_parcells.(names{i}) = pixels2rois(dFoF.(names{i}), parcells);
end

save(fullfile(outputPath, 'final_dFoF.mat'), 'dFoF','R','C', '-v7.3');
save(fullfile(outputPath, 'final_dFoF_parcels.mat'), 'dFoF_parcells', '-v7.3');
end 

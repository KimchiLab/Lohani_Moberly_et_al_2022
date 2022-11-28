function [Allparcells, parcells] = getParcellsByLansAllansAtlas

parcellsMat = 'allanParcellationTiffs\parcells.mat';
if exist(parcellsMat, 'file')
   load(parcellsMat,  'Allparcells', 'parcells');
   return;
end
load('allregions.mat', 'allregions');
[~,TXT] = xlsread('allanParcellationTiffs\subregion_list.csv');
N_parcells = length(TXT)-1;
Allparcells = zeros(size(allregions));
indicators = zeros([size(allregions) N_parcells]);
names = cell(N_parcells, 1);
description= cell(N_parcells, 1);
files = dir('allanParcellationTiffs\*.tif');
for n = 1:length(files)
    a = imread(fullfile(files(n).folder, files(n).name));
    a = imresize(a, size(allregions));
    ind = find(strcmpi(TXT(:,1), files(n).name(1:end-4)));
    indicators(:,:,ind-1) = double(rgb2gray(a)>0);
    names{ind-1} = TXT{ind,1};
    description{ind-1} = TXT{ind,2};
    Allparcells = Allparcells + double(rgb2gray(a)>0)*n;
end
parcells.names = names;
parcells.description = description;
parcells.indicators = indicators;
save(parcellsMat,  'Allparcells', 'parcells');
end 
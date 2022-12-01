%% Figure 2:  Differential coupling of behavioral variables to cholinergic and neural activity across the neocortex
% Lohani et al., 2022
clear all
% LOAD RAW DATA MAT FILE
%load('C:\Users\FIGURE 2\Fig2RawData.mat') % set path
load('C:\Users\AHM\Documents\cosyne_figs_for_SL\FIGURE 2\Fig2RawData.mat');
%% Figure2a: Full model heatmaps

figure()
subplot(2,1,1)
title('ACh')
% grabs
grabs_full = grabs_xls_full(1:6,:);
grabs_full_mu = nanmean(grabs_full);
plot_heatmap(grabs_full_mu,0,0.3,parcells_new);%axis square
subplot(2,1,2)
title('RCaMP')
% rcamp
rcamp_full = rcamp_xls_full(1:6,:);
rcamp_full_mu = nanmean(rcamp_full);
plot_heatmap(rcamp_full_mu,0,0.3,parcells_new);%axis square
%% Figure 2b: Full model ordered by ant post position:

figure()
subplot(1,2,1)
grabs_full = grabs_xls_full(1:6,:);
grabs_full(:,[21:26]) = []; %46 parcels unordered
grabs_full_left = grabs_full(:,2:2:end);
grabs_full_right = grabs_full(:,1:2:end); %23 parcels unordered

% parcel names unordered:
unordered_parcel_names = parcells_new.names;
unordered_parcel_names([21:26 53:56]) = []; %46 parcels unordered
unordered_parcel_names_left = unordered_parcel_names(2:2:end);
unordered_parcel_names_right = unordered_parcel_names(1:2:end);

% order based on ant-post center of mass
[num, txt] = xlsread(parcel_position,1);
ordered_parcels = txt(2:end,1); %23 ordered parcels WITH SINGLE QOUTES
center_mass = num(1:end,1); %23 positions
ordered_parcels_no_quotes = strrep(ordered_parcels, char(39), '');
[~, idx] = ismember(ordered_parcels_no_quotes,unordered_parcel_names_right);
% IDX is the new ordering of the parcels

grabs_full_left_trans = grabs_full_left';
grabs_full_left_trans = grabs_full_left_trans(idx,:);
%get corr of average:
[rho_full_left, pval] = corr(center_mass, nanmean(grabs_full_left_trans,2), 'type', 'Spearman');

hold on
% add mean and best fit line
g_full_left_mean = mean(grabs_full_left_trans,2);
g_full_left_std = std(grabs_full_left_trans,[],2);
g_full_left_sem = g_full_left_std./(sqrt(6));
errorbar(1:23,g_full_left_mean,g_full_left_sem,'.k','CapSize',0,'LineWidth',0.3,'MarkerSize',4)

coeffs = polyfit(1:23, g_full_left_mean, 1);
fittedX = linspace(min(1:23), max(1:23), 100);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 0.5);

ylim([-0.05 0.35])
text(2, 0.3,strcat('rho: ',num2str(rho_full_left,'%4.2f'),'pval=',num2str(pval,'%6.5f')))
title('ACh')
box off
set(gca,'TickDir','out')

subplot(1,2,2)
%** REPEAT FOR RCAMP PLOT
rcamp_xls_full = xlsread('C:\Users\AHM\Documents\cosyne_figs_for_SL\grabs_rcamp_fullAndShuffledVars.xlsx',5);

rcamp_full = rcamp_xls_full(1:6,:);
rcamp_full(:,[21:26]) = []; %46 parcels unordered
rcamp_full_left = rcamp_full(:,2:2:end);
rcamp_full_right = rcamp_full(:,1:2:end); %23 parcels unordered

% parcel names unordered:
unordered_parcel_names = parcells_new.names;
unordered_parcel_names([21:26 53:56]) = []; %46 parcels unordered
unordered_parcel_names_left = unordered_parcel_names(2:2:end);
unordered_parcel_names_right = unordered_parcel_names(1:2:end);

[num, txt] = xlsread(parcel_position,1);
ordered_parcels = txt(2:end,1); %23 ordered parcels WITH SINGLE QOUTES
center_mass = num(1:end,1); %23 positions
ordered_parcels_no_quotes = strrep(ordered_parcels, char(39), '');
[tf, idx] = ismember(ordered_parcels_no_quotes,unordered_parcel_names_right);
% IDX is the new ordering of the parcels

rcamp_full_left_trans = rcamp_full_left';
rcamp_full_left_trans = rcamp_full_left_trans(idx,:);
%get corr of average:
[rho_full_left, pval] = corr(center_mass, nanmean(rcamp_full_left_trans,2), 'type', 'Spearman');

hold on
% add mean and best fit line
r_full_left_mean = mean(rcamp_full_left_trans,2);
r_full_left_std = std(rcamp_full_left_trans,[],2);
r_full_left_sem = r_full_left_std./(sqrt(6));
errorbar(1:23,r_full_left_mean,r_full_left_sem,'.k','CapSize',0,'LineWidth',0.3,'MarkerSize',4)

coeffs = polyfit(1:23, r_full_left_mean, 1);
fittedX = linspace(min(1:23), max(1:23), 100);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 0.5);

ylim([-0.05 0.35])
text(2, 0.3,strcat('rho: ',num2str(rho_full_left,'%4.2f'),'pval=',num2str(pval,'%6.5f')))
title('rcamp')
box off
set(gca,'TickDir','out')

%% Figure 2c: Full and single variable mods

% full
grabs_full = grabs_xls_full(:,2:2:end);
grabs_full(grabs_full<0) = 0;
grabs_full_mu = nanmean(grabs_full,2);
% face 
grabs_fc = grabs_xls_face(:,2:2:end);
grabs_fc(grabs_fc<0) = 0;
grabs_fc_mu = nanmean(grabs_fc,2);
% pupil
grabs_pl = grabs_xls_pupil(:,2:2:end);
grabs_pl(grabs_pl<0) = 0;
grabs_pl_mu = nanmean(grabs_pl,2);
% wheel
grabs_wh = grabs_xls_wheel(:,2:2:end);
grabs_wh(grabs_wh<0) = 0;
grabs_wh_mu = nanmean(grabs_wh,2);

figure();
subplot(1,2,1);
hold on; title('Ach3.0')
scatter([0 0 0 0 0 0],grabs_full_mu,'.k');%plot(1,nanmean(grabs_full_mu),'r+')
full_modMu = nanmean(grabs_full_mu);
fullStd = nanstd(grabs_full_mu);
fullSEM = fullStd./sqrt(6);
bar(0,full_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(0,full_modMu,fullSEM,fullSEM,'r');

scatter([1 1 1 1 1 1],grabs_fc_mu,'.k');%plot(1,nanmean(grabs_full_mu),'r+')
fc_modMu = nanmean(grabs_fc_mu);
fcStd = nanstd(grabs_fc_mu);
fcSEM = fcStd./sqrt(6);
bar(1,fc_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(1,fc_modMu,fcSEM,fcSEM,'r');

scatter([2 2 2 2 2 2],grabs_pl_mu,'.k');%plot(1,nanmean(grabs_full_mu),'r+')
pl_modMu = nanmean(grabs_pl_mu);
plStd = nanstd(grabs_pl_mu);
plSEM = plStd./sqrt(6);
bar(2,pl_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(2,pl_modMu,plSEM,plSEM,'r');

scatter([3 3 3 3 3 3],grabs_wh_mu,'.k');%plot(1,nanmean(grabs_full_mu),'r+')
wh_modMu = nanmean(grabs_wh_mu);
whStd = nanstd(grabs_wh_mu);
whSEM = whStd./sqrt(6);
bar(3,wh_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(3,wh_modMu,whSEM,whSEM,'r');

ylim([-0.05 0.3]);
xticks([ 0 1 2 3]); xticklabels({'full' 'face' 'pupil' 'wheel'})

% repeat for RCAMP data
rcamp_xls_full = xlsread('C:\Users\AHM\Documents\cosyne_figs_for_SL\grabs_rcamp_fullAndShuffledVars.xlsx',5);
rcamp_xls_face = xlsread('C:\Users\AHM\Documents\cosyne_figs_for_SL\7_19_singleVar_grab_and_rcamp.xlsx',4);
rcamp_xls_pupil = xlsread('C:\Users\AHM\Documents\cosyne_figs_for_SL\7_19_singleVar_grab_and_rcamp.xlsx',5);
rcamp_xls_wheel = xlsread('C:\Users\AHM\Documents\cosyne_figs_for_SL\7_19_singleVar_grab_and_rcamp.xlsx',6);

% full
rcamp_full = rcamp_xls_full(:,2:2:end);
rcamp_full(rcamp_full<0) = 0;
rcamp_full_mu = nanmean(rcamp_full,2);
% face 
rcamp_fc = rcamp_xls_face(:,2:2:end);
rcamp_fc(rcamp_fc<0) = 0;
rcamp_fc_mu = nanmean(rcamp_fc,2);
% pupil
rcamp_pl = rcamp_xls_pupil(:,2:2:end);
rcamp_pl(rcamp_pl<0) = 0;
rcamp_pl_mu = nanmean(rcamp_pl,2);
% wheel
rcamp_wh = rcamp_xls_wheel(:,2:2:end);
rcamp_wh(rcamp_wh<0) = 0;
rcamp_wh_mu = nanmean(rcamp_wh,2);

subplot(1,2,2);
hold on; title('RCaMP')
scatter([0 0 0 0 0 0],rcamp_full_mu,'.k');
full_modMu = nanmean(rcamp_full_mu);
fullStd = nanstd(rcamp_full_mu);
fullSEM = fullStd./sqrt(6);
bar(0,full_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(0,full_modMu,fullSEM,fullSEM,'r');

scatter([1 1 1 1 1 1],rcamp_fc_mu,'.k');
fc_modMu = nanmean(rcamp_fc_mu);
fcStd = nanstd(rcamp_fc_mu);
fcSEM = fcStd./sqrt(6);
bar(1,fc_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(1,fc_modMu,fcSEM,fcSEM,'r');

scatter([2 2 2 2 2 2],rcamp_pl_mu,'.k');%plot(1,nanmean(grabs_full_mu),'r+')
pl_modMu = nanmean(rcamp_pl_mu);
plStd = nanstd(rcamp_pl_mu);
plSEM = plStd./sqrt(6);
bar(2,pl_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(2,pl_modMu,plSEM,plSEM,'r');

scatter([3 3 3 3 3 3],rcamp_wh_mu,'.k');%plot(1,nanmean(grabs_full_mu),'r+')
wh_modMu = nanmean(rcamp_wh_mu);
whStd = nanstd(rcamp_wh_mu);
whSEM = whStd./sqrt(6);
bar(3,wh_modMu,'FaceColor','none','EdgeColor','k');
er = errorbar(3,wh_modMu,whSEM,whSEM,'r');
ylim([-0.05 0.3])
xticks([ 0 1 2 3]); xticklabels({'full' 'face' 'pupil' 'wheel'})
%******************************************************************************************
% plots for deltas:
% plot full model and single variable model R2ds

grabs_noPupil = grabs_xls_noPupil(:,2:2:end);
grabs_noPupil(grabs_noPupil<0) = 0;
grabs_noPupil_mu = nanmean(grabs_noPupil,2);

grabs_noFace = grabs_xls_noFace(:,2:2:end);
grabs_noFace(grabs_noFace<0) = 0;
grabs_noFace_mu = nanmean(grabs_noFace,2);

grabs_noWheel = grabs_xls_noWheel(:,2:2:end);
grabs_noWheel(grabs_noWheel<0) = 0;
grabs_noWheel_mu = nanmean(grabs_noWheel,2);

noFace = grabs_full_mu-grabs_noFace_mu;
noFaceMu = nanmean(noFace);
noFaceStd = nanstd(noFace);
noFaceSEM = noFaceStd./sqrt(6); 

noPupil = grabs_full_mu-grabs_noPupil_mu;
noPupilMu = nanmean(noPupil);
noPupilStd = nanstd(noPupil);
noPupilSEM = noPupilStd./sqrt(6); 

noWheel = grabs_full_mu-grabs_noWheel_mu;
noWheelMu = nanmean(noWheel);
noWheelStd = nanstd(noWheel);
noWheelSEM = noWheelStd./sqrt(6); 

figure();
subplot(1,2,1);
hold on; title('Ach3.0')

scatter([1 1 1 1 1 1],grabs_full_mu-grabs_noFace_mu,'.k');%plot(2,nanmean(grabs_full_mu-grabs_noPupil_mu),'r+')
bar(1,noFaceMu,'FaceColor','none','EdgeColor','k');
er = errorbar(1,noFaceMu,noFaceSEM,noFaceSEM,'r'); 

scatter([2 2 2 2 2 2],grabs_full_mu-grabs_noPupil_mu,'.k');%plot(3,nanmean(grabs_full_mu-grabs_noFace_mu),'r+')
bar(2,noPupilMu,'FaceColor','none','EdgeColor','k');
er = errorbar(2,noPupilMu,noPupilSEM,noPupilSEM,'r'); 

scatter([3 3 3 3 3 3],grabs_full_mu-grabs_noWheel_mu,'.k');%plot(4,nanmean(grabs_full_mu-grabs_noWheel_mu),'r+')
bar(3,noWheelMu,'FaceColor','none','EdgeColor','k');
er = errorbar(3,noWheelMu,noWheelSEM,noWheelSEM,'r'); 
set(gca,'tickdir','out')
xticks([1 2 3])
xticklabels({'delta face' 'delta pupil' 'delta wheel'})
ylim([-0.05 0.3])

%*************************************************************************
% RCAMP PLOTS
% ***************************************************************
% plot average single variable model R2ds
rcamp_noPupil = rcamp_xls_noPupil(:,2:2:end);
rcamp_noPupil(rcamp_noPupil<0) = 0;
rcamp_noPupil_mu = nanmean(rcamp_noPupil,2);

rcamp_noFace = rcamp_xls_noFace(:,2:2:end);
rcamp_noFace(rcamp_noFace<0) = 0;
rcamp_noFace_mu = nanmean(rcamp_noFace,2);

rcamp_noWheel = rcamp_xls_noWheel(:,2:2:end);
rcamp_noWheel(rcamp_noWheel<0) = 0;
rcamp_noWheel_mu = nanmean(rcamp_noWheel,2);

noFace = rcamp_full_mu-rcamp_noFace_mu;
noFaceMu = nanmean(noFace);
noFaceStd = nanstd(noFace);
noFaceSEM = noFaceStd./sqrt(6); 

noPupil = rcamp_full_mu-rcamp_noPupil_mu;
noPupilMu = nanmean(noPupil);
noPupilStd = nanstd(noPupil);
noPupilSEM = noPupilStd./sqrt(6); 

noWheel = rcamp_full_mu-rcamp_noWheel_mu;
noWheelMu = nanmean(noWheel);
noWheelStd = nanstd(noWheel);
noWheelSEM = noWheelStd./sqrt(6); 

subplot(1,2,2);
hold on; title('jRCaMP1b')
full_modMu = nanmean(rcamp_full_mu);
fullStd = nanstd(rcamp_full_mu);
fullSEM = fullStd./sqrt(6);

scatter([1 1 1 1 1 1],rcamp_full_mu-rcamp_noFace_mu,'.k');%plot(2,nanmean(rcamp_full_mu-rcamp_noPupil_mu),'r+')
bar(1,noFaceMu,'FaceColor','none','EdgeColor','k');
er = errorbar(1,noFaceMu,noFaceSEM,noFaceSEM,'r'); 

scatter([2 2 2 2 2 2],rcamp_full_mu-rcamp_noPupil_mu,'.k');%plot(3,nanmean(rcamp_full_mu-rcamp_noFace_mu),'r+')
bar(2,noPupilMu,'FaceColor','none','EdgeColor','k');
er = errorbar(2,noPupilMu,noPupilSEM,noPupilSEM,'r'); 

scatter([3 3 3 3 3 3],rcamp_full_mu-rcamp_noWheel_mu,'.k');%plot(4,nanmean(rcamp_full_mu-rcamp_noWheel_mu),'r+')
bar(3,noWheelMu,'FaceColor','none','EdgeColor','k');
er = errorbar(3,noWheelMu,noWheelSEM,noWheelSEM,'r'); 
set(gca,'tickdir','out')
xticks([1 2 3])
xticklabels({'delta face' 'delta pupil' 'delta wheel'})
ylim([-0.05 0.3])

%%
function[] = plot_heatmap(correlation_vector, colorbar_low, colorbar_high,parcells_new)


for i = 1:1:56
    hold on;
    parcells=parcells_new.indicators;
    E=parcells(:,:,i); % overlay V1 parcell boundaries
    [B,~] = bwboundaries(E);
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 0.25)
    end
end

cm2 = vals2colormap(correlation_vector,'parula',[colorbar_low colorbar_high]);

for i = 1:1:56
    E=parcells(:,:,i);
    [B,~] = bwboundaries(E);
    for k = 1:length(B)
        boundary = B{k};
        
% remove the NAN parcels (21-26 and 53-56)

        if i == 21 || i == 22 || i == 23 || i == 24 || i == 25 || i == 26 || ...
        i == 53 || i == 54 || i == 55 || i == 56
        fill(boundary(:,2), boundary(:,1), 'w')
        else
        fill(boundary(:,2), boundary(:,1), cm2(i,:))
        end
    end
end

set(gca, 'box','on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
c = colorbar;
v = linspace(colorbar_low, colorbar_high, length(get(c,'TickLabels')));
set(c,'TickLabels', v);

end
%%
function rgb = vals2colormap(vals, colormap, crange)
% Take in a vector of N values and return and return a Nx3 matrix of RGB
% values associated with a given colormap
%
% rgb = AFQ_vals2colormap(vals, [colormap = 'jet'], [crange])
%
% Inputs:
% vals     = A vector of values to map to a colormap or a cell array of
%            vectors of values
% colormap = A matlab colormap. Examples: colormap = 'autumn';
%            colormap = 'jet'; colormap = 'hot';
% crange   = The values to map to the minimum and maximum of the colormap.
%            Defualts to the full range of values in vals.
%
% Outputs:
% rgb      = Nx3 matrix of rgb values mapping each value in vals to the
%            corresponding rgb colors.  If vals is a cell array then rgb
%            will be a cell array of the same length
%
% Example:
% vals = rand(1,100);
% rgb = AFQ_vals2colormap(vals, 'hot');
%
% Copyright Jason D. Yeatman, June 2012

if ~exist('colormap','var') || isempty(colormap)
    colormap = 'jet';
end

%
if ~iscell(vals)
    if ~exist('crange','var') || isempty(crange)
        crange = [min(vals) max(vals)];
    end
    % Generate the colormap
    cmap = eval([colormap '(256)']);
    % Normalize the values to be between 1 and 256
    vals(vals < crange(1)) = crange(1);
    vals(vals > crange(2)) = crange(2);
    valsN = round(((vals - crange(1)) ./ diff(crange)) .* 255)+1;
    % Convert any nans to ones
    valsN(isnan(valsN)) = 1;
    % Convert the normalized values to the RGB values of the colormap
    rgb = cmap(valsN, :);
elseif iscell(vals)
    if ~exist('crange','var') || isempty(crange)
        crange = [min(vertcat(vals{:})) max(vertcat(vals{:}))];
    end
    % Generate the colormap
    cmap = eval([colormap '(256)']);
    for ii = 1:length(vals)
        % Normalize the values to be between 1 and 256 for cell ii
        valsN = vals{ii};
        valsN(valsN < crange(1)) = crange(1);
        valsN(valsN > crange(2)) = crange(2);
        valsN = round(((valsN - crange(1)) ./ diff(crange)) .* 255)+1;
        % Convert any nans to ones
        valsN(isnan(valsN)) = 1;
        % Convert the normalized values to the RGB values of the colormap
        rgb{ii} = cmap(valsN, :);
    end
end
return
end
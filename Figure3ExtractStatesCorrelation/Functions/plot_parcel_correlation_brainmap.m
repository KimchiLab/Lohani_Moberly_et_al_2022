function[h]=plot_parcel_correlation_brainmap(correlation_vector, colormap,colorrange,CombinedParcellIdx,V1Idx,M2Idx,S1bIdx,pval)
%This function makes brain color maps from values in the correlation_vector
h=figure;
load('parcells_updated121519.mat'); parcells=parcells_new.indicators;
%convert values in correlation_vector to color map index
cm2 = vals2colormap(correlation_vector,colormap,colorrange);
%loop through each parcell and fill color between parcell boundaries 
for i = 1:length(CombinedParcellIdx)
    hold on; 
    E=parcells(:,:,CombinedParcellIdx(i));
    [B,~] = bwboundaries(E);
    for k = 1:length(B)
        boundary = B{k};
        if exist('pval') && pval(i) >=0.05 %if stats p-values are available
         fill(boundary(:,2), boundary(:,1), [211/255 211/255 211/255],'EdgeColor','none')
        else 
        fill(boundary(:,2), boundary(:,1), cm2(i,:),'EdgeColor','none')
        end 
    end
end
set(gca, 'colormap', str2num(colormap));
c = colorbar;
v = linspace(colorrange(1), colorrange(2), length(get(c,'TickLabels')));
set(c,'TickLabels', v);

E=parcells(:,:,CombinedParcellIdx(V1Idx)); % overlay V1 parcell boundaries
[B,~] = bwboundaries(E);
plot(B{1,1}(:,2), B{1,1}(:,1), 'Color',[231/255 49/255 51/255], 'LineWidth', 2)
E=parcells(:,:,CombinedParcellIdx(S1bIdx)); % overlay S1b parcell boundaries
B = bwboundaries(E);
plot(B{1,1}(:,2), B{1,1}(:,1), 'Color',[255/255 140/255 25/255], 'LineWidth', 2)
E=parcells(:,:,CombinedParcellIdx(M2Idx)); % overlay M2 parcell boundaries
B = bwboundaries(E);
plot(B{1,1}(:,2), B{1,1}(:,1),'Color',[210/255 0/255 255/255], 'LineWidth', 2)

axis off; box off; 

end



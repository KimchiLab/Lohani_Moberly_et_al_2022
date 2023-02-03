function [h1]=plot_parcell_overlay(mov,R,C,plotSave,parcells)
if plotSave
    meanx=nanmean(mov,2);
    mov_trans=reshape(meanx,R,C);
    h1=figure; subplot(1,2,1);imshow(imadjust(mat2gray(mov_trans)));
    subplot(1,2,2);imshow(imadjust(mat2gray(mov_trans)));hold on; 
    for k=1:size(parcells,3)
        E=parcells(:,:,k);
        [B,~] = bwboundaries(E);hold on;
        for t = 1:length(B)
            boundary = B{t};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.1)
        end
    end
    hold off;
else
    h1=[];
end

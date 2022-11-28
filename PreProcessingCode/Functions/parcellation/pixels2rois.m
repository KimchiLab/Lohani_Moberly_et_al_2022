 function roi_time_series = pixels2rois(braincs_WB, parcells)
%% average over all active pixels in a parcell 
[~, N_time, N_trials]  = size(braincs_WB);

roi_time_series = zeros(length(parcells.names),  N_time, N_trials);

for pi = 1:length(parcells.names)    
    pI = find(parcells.indicators(:,:,pi)>0); 
    pI=pI(pI<=size(braincs_WB,1));
    currTimeSeries=braincs_WB(pI, :,:);
    indices=find(max(abs(currTimeSeries),[],2)>0);%remove pixels with no activity (this will also remove all saturated pixels)
    if (numel(indices)/numel(pI))>0.25 %get a mean across all pixels only if at least 25% of pixels in the parcell are active, otherwise throw out the entire parcell  
    currTimeSeries=currTimeSeries(indices,:,:); 
    roi_time_series(pi, :, :) = mean(currTimeSeries,1);   
    else
    roi_time_series(pi, :, :) = nan;  
    end 
    clear currTimeSeries indices 
end

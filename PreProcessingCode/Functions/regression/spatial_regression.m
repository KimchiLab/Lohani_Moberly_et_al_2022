function [output_sig, sig_blue, sig_uv] = spatial_regression(blue_sig, uvsig, brain_mask, patch_size, delay_length) 

% Regression of hemodynamics, exploiting spatial info.
% Input: 
%  blue_sig     - calcium dependent signal - pixels X time (not masked)
%  uv_sig       - calcium indipendendt - pixels over time (not masked)
%  brain_mask   - matrix of R*C = #pixels of zeors and ones indicating
%                 location of brain.
%  patchsize    - patch size for spatial info, preferably mor
%  delay_length - ind of time frame for which detrending/filtering ends
%
% Output: 
%   output_sig  - regressed signal, pixels X time (nans where brain_mask = 0
%   sig_blue    - estimated standard deviation of noise in blue - pixels X 1
%   sig_uv      - estimated standard deviation of noise in uv - pixels X 1
%     
% Written by:
%    Hadas Benisty, hadas.benisty@gmail.com
%    Boris Landa
%    April 6 2020
%
%------------------------------------------------------------------------

% this analysis requires zero mean signals so removing them and later
% adding them
mean_blue =  mean(blue_sig(:,delay_length+1:end),2);
mean_uv =  mean(uvsig(:,delay_length+1:end),2);

rmmean_blue = bsxfun(@minus, blue_sig, mean_blue);
rmmean_uv = bsxfun(@minus, uvsig, mean_uv);
patch_size = round((patch_size+1)/2);

[sig_blue, sig_uv, output_sig]  = regress_filter(rmmean_blue(:, delay_length+1:end), rmmean_uv(:, delay_length+1:end), brain_mask, patch_size);
output_sig = cat(2, nan(size(output_sig, 1), delay_length), output_sig); 
output_sig=output_sig+mean_blue; %adding the mean singal back , SL added this 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Detrending of time series of a single pixel using 2 methods
% Input: 
%   - insig - time series of a pixel
%   - method - 'FIR' (default) or 'TOPHAT'
%   - filtLen - filter length (default = 100 for FIR, 1000 for TOPHAT)
%   - filtcutoff - used for FIR, cuttoff frequency (default = 0.01)
%
%  Written by Hadas Benisty 11/11/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputsig = detrendPixel(insig, method, filtLen, filtcutoff)

if nargin == 1
    filtLen = 100;filtcutoff = 0.01;
    method = 'FIR';
end
switch method
    case 'TOPHAT'
        if ~exist('filtLen', 'var')
            filtLen = 1e3;
        end
        se = strel('line', filtLen, 0);
        outputsig=imtophat([2*insig(1)-fliplr(insig(1:filtLen)),insig],se);
        outputsig=outputsig(filtLen+1:end)';
    case 'FIR'        
        filtcoeff = fir1(filtLen,filtcutoff);
        baseline = filter(filtcoeff, 1, [2*insig(1)-fliplr(insig(1:filtLen)),insig]);    
        baseline = baseline(filtLen+1:end);
        outputsig=insig - baseline;
        outputsig=outputsig-min(outputsig);
end
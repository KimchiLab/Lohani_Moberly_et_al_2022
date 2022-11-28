function [sig_noise1, sig_noise2, data_filt] = regress_filter(blue_sig, uv_sig, brain_mask, patchsize)

% Regression of hempdynamics based on spatial info.
% Input: 
%  blue_sig    - calcium dependent signal - pixels X time (not masked)
%  uv_sig      - calcium indipendendt - pixels over time (not masked)
%  brain_mask  - matrix of R*C = #pixels of zeors and ones indicating
%               location of brain.
% patchsize    - half size of patch + 1 for spatial info
%
% Output: 
%   sig_noise1 - estimated standard deviation of noise in blue - pixels X 1
%   sig_noise1 - estimated standard deviation of noise in uv - pixels X 1
%   data_filt  - regressed signalj
%     
% Written by:
%    Hadas Benisty, hadas.benisty@gmail.com
%    Boris Landa
%    April 6 2020
%
%------------------------------------------------------------------------

disp(['Begining of filtering + regression']);

[R,C] = size(brain_mask);
[Y,X]=meshgrid(1:R,1:C);
X=X(:);Y=Y(:);
braininds = find(brain_mask);
tic;
[Np,Nt] = size(blue_sig);
sig_noise1 = zeros(Np,1);
sig_noise2 = zeros(Np,1);
data_filt = nan(Np, Nt);
%cluster = parcluster;
for pii=1:length(X)  
    if mod(pii,1000)==0
        disp(['Filtering + Regression: processed ' num2str(pii/length(X)) ' pixels']);
    end
    x = X(pii);y = Y(pii);
    nn = find(abs(X-x) < patchsize & abs(Y-y) < patchsize);
    
    nn0 = find(abs(X-x) ==0 & abs(Y-y) ==0);
    if ismember(nn0 , braininds)  %SL replaced all continue statememnts with if stamements
        
        nn = intersect(nn, braininds);
        if ~isempty(nn)
            
            midpixel = find(nn==nn0);
            
            y1 = blue_sig(nn, :);
            y2 = uv_sig(nn, :);
            
            naninds=find(isnan(y1(:,1))); %SL added this part
            naninds1=find(isnan(y2(:,1))); 
            y1(naninds,:)=0; y2(naninds,:)=0; %SL added this part
            y1(naninds1,:)=0; y2(naninds1,:)=0; %SL added this part
            
            if ~(all(y1(:)==0)||  all(y2(:)==0)) % skip if values are all empty , SL added this 
               
            % evaluating noise sigma
            % to do - evalute for more experiments and see what's consistent
            %     sig_noise1(pii) = estimateNoiseCov(y1);
            %     sig_noise2(pii) = estimateNoiseCov(y2);
            %     sig_noise1=1.9e3;
            %     sig_noise2=1.5e3;
            d = size(y1,1);
            y = cat(1, y1, y2);
            sigy = y*y.'/size(y,2);
            sigy1 = sigy(1:d,1:d);
            sigy12 = sigy(1:d,d+1:end);
            sigy2 = sigy(1+d:end,1+d:end);
            
            sig_noise1(pii) = median(svd(sigy1));
            sig_noise2(pii) = median(svd(sigy2));
            
            sigx = sigy1-sig_noise1(pii)*eye(d) - sigy12*pinv((sigy2-sig_noise2(pii)*eye(d)))*sigy12';
            % to do - actually we can save time and estimate xhat only for the pixel in the
            % middle of the patch
            xhat = [sigx zeros(d)]*pinv(sigy)*cat(1, blue_sig(nn, :), uv_sig(nn, :));
            data_filt(pii,:) = xhat(midpixel,:);
            end 
        end
    end
end
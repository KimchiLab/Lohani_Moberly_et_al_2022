VisualPixel_2d=[82,174]; 
VisPix=sub2ind([R, C],VisualPixel_2d(1),VisualPixel_2d(2));
patch_size=14; delay_length=0; 
mean_blue =  mean(dFoF.blue(:,delay_length+1:end),2);
mean_uv =  mean(dFoF.uv(:,delay_length+1:end),2);

rmmean_blue = bsxfun(@minus, dFoF.blue, mean_blue);
rmmean_uv = bsxfun(@minus, dFoF.blue, mean_uv);
patchsize = round((patch_size+1)/2);

load('brainMask.mat','mask');
maskinds = find(mask);
naninds=find(isnan(dFoF.blue(:,1)));
for i=1:length(naninds)
    [naninds_sub(i,1),naninds_sub(i,2)]=ind2sub([R,C],naninds(i));
    mask(naninds_sub(i,1),naninds_sub(i,2))=0;
end
brain_mask=mask; 

blue_sig=rmmean_blue; uv_sig=rmmean_uv; 
[R,C] = size(brain_mask);
[Y,X]=meshgrid(1:R,1:C);
X=X(:);Y=Y(:);
braininds = find(brain_mask);
[Np,Nt] = size(blue_sig);
sig_noise1 = zeros(Np,1);
sig_noise2 = zeros(Np,1);

    x = VisualPixel_2d(1);y = VisualPixel_2d(2);
    nn = find(abs(X-x) < patchsize & abs(Y-y) < patchsize);
    
    nn0 = find(abs(X-x) ==0 & abs(Y-y) ==0);
    nn = intersect(nn, braininds);
    midpixel = find(nn==nn0);
     
    y1 = blue_sig(nn, :);
    y2 = uv_sig(nn, :);
    
    d = size(y1,1);
    y = cat(1, y1, y2);
    sigy = y*y.'/size(y,2);
    sigy1 = sigy(1:d,1:d);
    sigy12 = sigy(1:d,d+1:end);
    sigy2 = sigy(1+d:end,1+d:end);
    
     sig_noise1(pii) = median(svd(sigy1));
     sig_noise2(pii) = median(svd(sigy2));
    
    sigx = sigy1-sig_noise1(pii)*eye(d) - sigy12*pinv((sigy2-sig_noise2(pii)*eye(d)))*sigy12';
  xhat = [sigx zeros(d)]*pinv(sigy)*cat(1, blue_sig(nn, :), uv_sig(nn, :));
    data_filt(pii,:) = xhat(midpixel,:);
    
    
    H=[sigx zeros(d)]*pinv(sigy); 
    H_fin=H(midpixel,1:d);
    H_fin1=H(midpixel,d+1:end);
    final_data=nan(R*C,1); 
    final_data(nn)=H_fin; 
    figure;imagesc(reshape(final_data,R,C));xlim([VisualPixel_2d(2)-patchsize  VisualPixel_2d(2)+patchsize]);  ylim([VisualPixel_2d(1)-patchsize  VisualPixel_2d(1)+patchsize]);    
  
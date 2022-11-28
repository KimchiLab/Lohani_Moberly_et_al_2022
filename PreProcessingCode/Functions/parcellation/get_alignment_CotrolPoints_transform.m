function [tform,R,C] = get_alignment_CotrolPoints_transform(insig,template,R,C,movingPoints,fixedPoints)
%does control point matching by selecting points on moving image (insig) a
% and template(parcellated template) . Need to select at least three points
toalignframe=reshape(insig(:,1),R,C); 
toalignframe=imadjust(mat2gray(toalignframe),[0.45 0.95],[]);
disp('Please select control points on the brain and close GUI when done');
[movingPoints, fixedPoints]=cpselect(toalignframe,template,movingPoints,fixedPoints,'Wait',true);%call builit in GUI to select points 
tform = fitgeotrans(movingPoints,fixedPoints,'Similarity');%similatiry based transformation
Jregistered = imwarp(toalignframe,tform,'OutputView',imref2d(size(template)));
figure;imshowpair(template,Jregistered,'diff')
end





function [tform,R,C,movingPoints,fixedPoints] = get_alignment_CotrolPoints_transform_prompt(insig,template,R,C,movingPoints,fixedPoints)
%does control point matching by selecting points on moving image (insig) a
% and template(parcellated template) . Need to select at least three points
str='Y';
brightIdx=[];
%while user is not satisfied with registration, keep asking if user want to redo registration
while strcmp(str,'Y')
    toalignframe=reshape(insig(:,1),R,C);
    if isempty(brightIdx)
        toalignframe=imadjust(mat2gray(toalignframe));
    else
        toalignframe=imadjust(mat2gray(toalignframe),brightIdx,[]);
    end
    disp('Please select control points on the brain and close GUI when done');
    [movingPoints, fixedPoints]=cpselect(toalignframe,template,movingPoints,fixedPoints,'Wait',true);%call builit in GUI to select points
    tform = fitgeotrans(movingPoints,fixedPoints,'Similarity');%similatiry based transformation
    Jregistered = imwarp(toalignframe,tform,'OutputView',imref2d(size(template)));
    figure;imshowpair(template,Jregistered,'diff')
    
    % check the output of registration amd allow interactive use if not good
    prompt = 'Do you want to redo registration? Y/N [Y]: ';  
     str = input(prompt,'s');
    if isempty(str)|| strcmp(str,'N')
        disp('Continuing on without additional modification to registration')
        break;        
    end
    prompt1='Adjust moving image brightness by inputing the idx between 0 and 1 as something like [0.3 0.7]:  ';   
    brightIdx=input(prompt1); 
end





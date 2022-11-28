function [outMov]=readCxdSingleFile(tiffsPath,resizeRatio)
tiffs=dir(fullfile(tiffsPath, '/*.cxd'));
filepath=fullfile(tiffsPath,tiffs.name);
data=bfopen(filepath); 
tmpMov=data{1,1}(:,1);
tmp2Mov=cell2mat(tmpMov');
T=size(tmpMov,1);[R,C]=size(tmpMov{1,1}); 
outMov=reshape(tmp2Mov,R,C,T);

if resizeRatio ~= 1
        outMov = imresize(outMov, resizeRatio);
end
end 
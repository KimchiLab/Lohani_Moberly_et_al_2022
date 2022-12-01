
function[aveData]=averageStateData(data_final,numParcells) 
aveData=nan(numParcells,size(data_final,1));
for jj=1:size(data_final,1) %for each animal 
    %combine across sessions
    currData=cell(1,size(data_final,2)); 
    for kk=1:size(data_final,2)%for each session 
        if ~isempty(data_final{jj,kk})
            currData{kk}=cell2mat(data_final{jj,kk});
        end
    end
    %combine data across epcohs and sessions and then take an average 
        currData1=cell2mat(currData); 
        if ~isempty(currData1)
            aveData(:,jj)=nanmean(currData1,2)*100;%multiply by 100 to get percentage 
        end
end
function[state1corrMat, state1corrMat_Shuffled,state2corrMat,state2corrMat_Shuffled,state3corrMat,state3corrMat_Shuffled,state1Imaging,state2Imaging,state3Imaging,state1Imaging_Shuffled,state2Imaging_Shuffled,state3Imaging_Shuffled]...
    =statePearsonCorrelation(state1On,state1Off,state2On,state2Off,state3On,state3Off,imaging_time, dFoF_parcells,dFoF_parcells_Shuffled,names,signalsExtraction)
%% run Pearson's correlations on state1,state2, and state3 data
%% inputs
%state1On: onset times for state1 (such as locomotion)
%state1Off: offset times for state1
%state2On: onset times for state2 (such as face high quiescence)
%state2Off: offset times for state2
%state3On: onset times for state3 (such as face low quiescence)
%state3Off: offset times for state3
%imaging_time: imaging timestamps 
%dFoF_parcells: imaging data by parcells
%dFoF_parcells_Shuffled: shuffled imaging data by parcells
%names: names of all imaging colors
%signalsExtraction: 'blueuv' or 'RCaMP_AC'  
%% outputs
%state1corrMat/state2corrMat/state3corrMat: correlation matrix for each state as numparcellsXnumParcellsx number of state
%state1corrMat_Shuffled/state2corrMat_Shuffled/state3corrMat_Shuffled: correlation matrix for shuffled data
%state1Imaging/state2Imaging/state3Imaging: imaging data organized around states
%state1Imaging_Shuffled/state2Imaging_Shuffled/state3Imaging_Shuffled: shuffled imaging data organized around states 
%%
numParcells=size(dFoF_parcells.blue,1);
%for each color 
for i=1:length(names)
    state1corrMat.(names{i})=zeros(numParcells,numParcells,length(state1On));  state1corrMat_Shuffled.(names{i})=zeros(numParcells,numParcells,length(state1On));
    state2corrMat.(names{i})=zeros(numParcells,numParcells,length(state2On));  state2corrMat_Shuffled.(names{i})=zeros(numParcells,numParcells,length(state2On));
    state3corrMat.(names{i})=zeros(numParcells,numParcells,length(state3On));  state3corrMat_Shuffled.(names{i})=zeros(numParcells,numParcells,length(state3On));
    
    %for each state, calculate parcell-wise Pearson's correlation for each region pair 
    for r=1:length(state1On)
        state1Imaging.(names{i}){r}=dFoF_parcells.(names{i})(:,(imaging_time>state1On(r) & imaging_time<state1Off(r)));
        state1corrMat.(names{i})(:,:,r)=corr(state1Imaging.(names{i}){r}');
        state1corrMat.(names{i})(:,:,r)=state1corrMat.(names{i})(:,:,r)-diag(diag(state1corrMat.(names{i})(:,:,r))); % make diagonal elements 0 becaue it's autocorrelation
        
        state2Imaging.(names{i}){r}=dFoF_parcells.(names{i})(:,(imaging_time>state2On(r) & imaging_time<state2Off(r)));
        state2corrMat.(names{i})(:,:,r)=corr(state2Imaging.(names{i}){r}');
        state2corrMat.(names{i})(:,:,r)=state2corrMat.(names{i})(:,:,r)-diag(diag(state2corrMat.(names{i})(:,:,r))); % make diagonal elements 0 becaue it's autocorrelation
        
        state3Imaging.(names{i}){r}=dFoF_parcells.(names{i})(:,(imaging_time>state3On(r) & imaging_time<state3Off(r)));
        state3corrMat.(names{i})(:,:,r)=corr(state3Imaging.(names{i}){r}');
        state3corrMat.(names{i})(:,:,r)=state3corrMat.(names{i})(:,:,r)-diag(diag(state3corrMat.(names{i})(:,:,r))); % make diagonal elements 0 becaue it's autocorrelation
        
        %do correlations for each surrogate/shuffled dataset and take a mean across datasets
        state1Imaging_Shuffled.(names{i}){r}=dFoF_parcells_Shuffled.(names{i})(:,(imaging_time>state1On(r) & imaging_time<state1Off(r)),:);
        tmp =zeros(numParcells,numParcells,size(state1Imaging_Shuffled.(names{i}){r},3));
        for numShuff=1:size(state1Imaging_Shuffled.(names{i}){r},3)
            tmp(:,:,numShuff) =corr(state1Imaging_Shuffled.(names{i}){r}(:,:,numShuff)');
            tmp(:,:,numShuff)=tmp(:,:,numShuff)-diag(diag(tmp(:,:,numShuff))); % make diagonal elements 0 becaue it's autocorrelation
        end
        state1corrMat_Shuffled.(names{i})(:,:,r)=mean(tmp,3);
        
        state2Imaging_Shuffled.(names{i}){r}=dFoF_parcells_Shuffled.(names{i})(:,(imaging_time>state2On(r) & imaging_time<state2Off(r)),:);
        tmp =zeros(numParcells,numParcells,size(state2Imaging_Shuffled.(names{i}){r},3));
        for numShuff=1:size(state2Imaging_Shuffled.(names{i}){r},3)
            tmp(:,:,numShuff) =corr(state2Imaging_Shuffled.(names{i}){r}(:,:,numShuff)');
            tmp(:,:,numShuff)=tmp(:,:,numShuff)-diag(diag(tmp(:,:,numShuff))); % make diagonal elements 0 becaue it's autocorrelation
        end
        state2corrMat_Shuffled.(names{i})(:,:,r)=mean(tmp,3);
        
        state3Imaging_Shuffled.(names{i}){r}=dFoF_parcells_Shuffled.(names{i})(:,(imaging_time>state3On(r) & imaging_time<state3Off(r)),:);
        tmp =zeros(numParcells,numParcells,size(state3Imaging_Shuffled.(names{i}){r},3));
        for numShuff=1:size(state3Imaging_Shuffled.(names{i}){r},3)
            tmp(:,:,numShuff) =corr(state3Imaging_Shuffled.(names{i}){r}(:,:,numShuff)');
            tmp(:,:,numShuff)=tmp(:,:,numShuff)-diag(diag(tmp(:,:,numShuff))); % make diagonal elements 0 becaue it's autocorrelation
        end
        state3corrMat_Shuffled.(names{i})(:,:,r)=mean(tmp,3);
    end
end
%get correlation between blue and green channcels if both colors exist
if strcmp(signalsExtraction,'RCaMP_AC')
    state1corrMat.BG=zeros(numParcells,numParcells,length(state1On));  state1corrMat_Shuffled.BG=zeros(numParcells,numParcells,length(state1On));
    state2corrMat.BG=zeros(numParcells,numParcells,length(state2On));  state2corrMat_Shuffled.BG=zeros(numParcells,numParcells,length(state2On));
    state3corrMat.BG=zeros(numParcells,numParcells,length(state3On));  state3corrMat_Shuffled.BG=zeros(numParcells,numParcells,length(state3On));
    
    for r=1:length(state1On)
        state1corrMat.BG(:,:,r)=corr(state1Imaging.blue{r}',state1Imaging.green{r}');
        state2corrMat.BG(:,:,r)=corr(state2Imaging.blue{r}',state2Imaging.green{r}');
        state3corrMat.BG(:,:,r)=corr(state3Imaging.blue{r}',state3Imaging.green{r}');
        
        %do correlations for each surrogate/shuffled dataset and take a mean across datasets
        tmp =zeros(numParcells,numParcells,size(state1Imaging_Shuffled.blue{r},3));
        for numShuff=1:size(state1Imaging_Shuffled.blue{r},3)
            tmp(:,:,numShuff) =corr(state1Imaging_Shuffled.blue{r}(:,:,numShuff)',state1Imaging_Shuffled.green{r}(:,:,numShuff)');
        end
        state1corrMat_Shuffled.BG(:,:,r)=mean(tmp,3);
        
        tmp =zeros(numParcells,numParcells,size(state2Imaging_Shuffled.blue{r},3));
        for numShuff=1:size(state2Imaging_Shuffled.blue{r},3)
            tmp(:,:,numShuff) =corr(state2Imaging_Shuffled.blue{r}(:,:,numShuff)',state2Imaging_Shuffled.green{r}(:,:,numShuff)');
        end
        state2corrMat_Shuffled.BG(:,:,r)=mean(tmp,3);
        
        tmp =zeros(numParcells,numParcells,size(state3Imaging_Shuffled.blue{r},3));
        for numShuff=1:size(state3Imaging_Shuffled.blue{r},3)
            tmp(:,:,numShuff) =corr(state3Imaging_Shuffled.blue{r}(:,:,numShuff)',state3Imaging_Shuffled.green{r}(:,:,numShuff)');
        end
        state3corrMat_Shuffled.BG(:,:,r)=mean(tmp,3);
    end
end
end

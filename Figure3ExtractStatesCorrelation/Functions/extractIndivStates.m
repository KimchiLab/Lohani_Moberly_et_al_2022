function[state1Imaging,state2Imaging,state3Imaging]=extractIndivStates(state1On,state1Off,state2On,state2Off,state3On,state3Off,imaging_time, dFoF_parcells,names)
%% extract state1,state2, and state3 imaging data
%% inputs
%state1On: onset times for state1 (such as locomotion)
%state1Off: offset times for state1
%state2On: onset times for state2 (such as face high quiescence)
%state2Off: offset times for state2
%state3On: onset times for state3 (such as face low quiescence)
%state3Off: offset times for state3
%imaging_time: imaging timestamps
%dFoF_parcells: imaging data by parcells
%names: names of all imaging colors
%% outputs
%state1Imaging/state2Imaging/state3Imaging: imaging data organized around states
%% extract data for each color
for i=1:length(names)
    state1Imaging.(names{i})=[]; state2Imaging.(names{i})=[]; state3Imaging.(names{i})=[];
    if ~isempty(state1On)
        for r=1:length(state1On)
            state1Imaging.(names{i}){r}=dFoF_parcells.(names{i})(:,(imaging_time>state1On(r) & imaging_time<state1Off(r)));
        end
    end
    
    if ~isempty(state2On)
        for r=1:length(state2On)
            state2Imaging.(names{i}){r}=dFoF_parcells.(names{i})(:,(imaging_time>state2On(r) & imaging_time<state2Off(r)));
        end
    end
    
    if ~isempty(state3On)
        for r=1:length(state3On)
            state3Imaging.(names{i}){r}=dFoF_parcells.(names{i})(:,(imaging_time>state3On(r) & imaging_time<state3Off(r)));
        end
    end
end
end

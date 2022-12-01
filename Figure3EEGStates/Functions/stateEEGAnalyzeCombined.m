function[state1EEG,state2EEG,state3EEG,state1EEG_env,state2EEG_env,state3EEG_env,state1_lowpower,state2_lowpower,state3_lowpower,state1_highpower,state2_highpower,state3_highpower]=stateEEGAnalyzeCombined...
    (state1On,state1Off,state2On,state2Off,state3On, state3Off,eeg_time, eeg_raw,eeg_env,Fs)
%% extract eeg bandpower during sustained states 
%% inputs
%state1On: onset times for state1 (such as locomotion)
%state1Off: offset times for state1
%state2On: onset times for state2 (such as face high)
%state2Off: offset times for state2
%state3On: onset times for state2 (such as face low)
%state3Off: offset times for state2
%eeg_time: eeg timestamps 
%eeg_raw: raw eeg
%eeg_env: eeg low frequency amplitude 
%Fs: EEG sampling frequency 

%%outputs
%state1EEG/state2EEG/state3EEG: EEG extracted during states 1, 2 and 3 
%state1EEG_env/state2EEG_env/state3EEG_env:  EEG low frequency envelope extracted during states 
%state1_lowpower/state2_lowpower/state3_lowpower: low frequency(1-10 Hz) band power 
%state1_highpower/state2_highpower/state3_highpower: high frequency(30-100 Hz) band power 

    state1EEG=cell(1,length(state1On));  
    state2EEG=cell(1,length(state1On)); 
    state3EEG=cell(1,length(state1On)); 
    state1EEG_env=nan(1,length(state1On));  
    state2EEG_env=nan(1,length(state1On)); 
    state3EEG_env=nan(1,length(state1On)); 
    state1_lowpower=nan(1,length(state1On)); 
    state2_lowpower=nan(1,length(state1On)); 
    state3_lowpower=nan(1,length(state1On)); 
    state1_highpower=nan(1,length(state1On)); 
    state2_highpower=nan(1,length(state1On)); 
    state3_highpower=nan(1,length(state1On)); 
    %for each state, extract raw EEG and calculate bandpower, also extract EEG low frequency amplitude, 
    for r=1:length(state1On)
        state1EEG{r}=eeg_raw(eeg_time>state1On(r) & eeg_time<state1Off(r));
        state2EEG{r}=eeg_raw(eeg_time>state2On(r) & eeg_time<state2Off(r));
        state3EEG{r}=eeg_raw(eeg_time>state3On(r) & eeg_time<state3Off(r));
        
        %calculate bandpower in high and low frequencies  
        state1_lowpower(r) = bandpower(state1EEG{r},Fs,[1 10]);
        state1_highpower(r) = bandpower(state1EEG{r},Fs,[30 100]);
        
        state2_lowpower(r) = bandpower(state2EEG{r},Fs,[1 10]);
        state2_highpower(r) = bandpower(state2EEG{r},Fs,[30 100]);
        
        state3_lowpower(r) = bandpower(state3EEG{r},Fs,[1 10]);
        state3_highpower(r) = bandpower(state3EEG{r},Fs,[30 100]);

        %extract low frequency amplitude and average the continuous amplitude signal 
        state1EEG_env(r)=mean(eeg_env(eeg_time>state1On(r) & eeg_time<state1Off(r)));
        state2EEG_env(r)=mean(eeg_env(eeg_time>state2On(r) & eeg_time<state2Off(r)));
        state3EEG_env(r)=mean(eeg_env(eeg_time>state3On(r) & eeg_time<state3Off(r)));
              
    end
end

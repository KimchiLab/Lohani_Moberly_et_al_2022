Code for figures 3,4, and 5 in https://www.biorxiv.org/content/10.1101/2020.12.09.418632v2 

First run AdjustFacePCNormalizeVideoVars on facemap proc outputs to flip the signs of PC1 if they are negatively skewed and normalize. 

Then run ParcellWiseCorrelation.m for doing either "pearson's non-lagged" correlations or "cross-correlations" on sustained locomotion, face high, and face low states. The sustained states are set to be of 5s duration (this was so that we could get enough frames for cross-correlations). If you are extracting state data for Pearson's correlation or for state amplitude, you can set it to < 5s. 
refer to Figures 4 and 5 in the paper

The function  CharacterizeSustainedStateAmplitude gets the mean df/f in sustained locomotion,face high,
and face low states and does statistics. 
refer to Figures3b and S5d in the paper. 

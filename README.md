MATLAB code for preprocessing dual-color mesoscopic images, data analysis, and generating figures for the paper Lohani & Moberly et al. (2022). Spatiotemporally heterogeneous coordination of cholinergic and neocortical activity. https://www.nature.com/articles/s41593-022-01202-6

1. Preprocessing Code: Preprocess raw tif (or cxd) video files, sync behavior, and extract uv corrected df/f
2. Figure2Regression: Cross-validated explained variance (R2) based on a ridge regression model fitted with FaceMap PCs, locomotion and pupil area
3. Figure3EEGStates: EEG time-frequency analysis for three sustained behavioral states (face low, face high, locomotion)
4. Figure3ExtractStatesCorrelation:Correlations in three sustained behavioral states between parcells or within parcell for two colors
5. FigureS3HemodynamicCorrection:Peak of negative flourescence dips (hemodynamic artifact) around a particular event
6. FigureS7TransientStates: Event triggered averages at state transitions

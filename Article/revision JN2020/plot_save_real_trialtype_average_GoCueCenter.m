% plot_save_real_trialtype_average_GoCueCenter Script JC 101/02/2020
% load('Ntrial_type.mat')

idx1 = trial.idx_correct_L
idx2 = trial.idx_correct_R
idx_str1 = 'corL'
idx_str2 = 'corR'

% VAR_DLC=MouthY_tr;
% VAR_str = 'MouthY'
% plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);
% % saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.emf']);

VAR_DLC=NoseY_tr;
VAR_str = 'NoseY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.emf']);

% VAR_DLC=NoseX_tr;
% VAR_str = 'NoseX'
% plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);
% % saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.emf']);


VAR_DLC=TongueY_tr;
VAR_str = 'TongueY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_' MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_' MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.emf']);

% VAR_DLC=TongueX_tr;
% VAR_str = 'TongueX'
% plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_' MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);

% VAR_DLC=WhiskerY_tr;
% VAR_str = 'WhiskerY'
% plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_' MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_' MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.emf']);

VAR_DLC=WhiskerX_tr;
VAR_str = 'WhiskerX'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_' MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_' MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.emf']);

%% omi vs omiopto trials selection
idx1 = trial.idx_NoLick
idx2 = trial.idx_NoLick_opto
idx_str1 = 'omi'
idx_str2 = 'omi-opto'

% VAR_DLC=MouthY_tr;
% VAR_str = 'MouthY'
% plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);

VAR_DLC=NoseY_tr;
VAR_str = 'NoseY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);

VAR_DLC=TongueY_tr;
VAR_str = 'TongueY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);

%% omi vs cor trials selection
idx1 = [trial.idx_correct_R trial.idx_correct_L]
idx2 = trial.idx_NoLick
idx_str1 = 'cor'
idx_str2 = 'omi'

% VAR_DLC=MouthY_tr;
% VAR_str = 'MouthY'
% plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);

VAR_DLC=NoseY_tr;
VAR_str = 'NoseY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);

VAR_DLC=TongueY_tr;
VAR_str = 'TongueY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);

%% cor vs opto trials selection
idx1 = [trial.idx_correct_L trial.idx_correct_R]
idx2 = trial.idx_all_opto
idx_str1 = 'cor-all'
idx_str2 = 'opto-all'

% VAR_DLC=MouthY_tr;
% VAR_str = 'MouthY'
% plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);


VAR_DLC=NoseY_tr;
VAR_str = 'NoseY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);


VAR_DLC=TongueY_tr;
VAR_str = 'TongueY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
saveas(gcf,[SAVEFIG_folder 'DLC_TrialAverage_'  MouseID '_' idx_str1 '_' idx_str2  '_' VAR_str '.png']);
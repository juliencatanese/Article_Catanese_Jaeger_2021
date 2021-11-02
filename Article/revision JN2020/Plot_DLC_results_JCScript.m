% Plot_DLC_results_JCScript
% Julien Catanese 10/07/2020

SAVEFIG_folder = 'DLC_results\'
mkdir(['.\'  SAVEFIG_folder ])

load([ MouseID '_DLCresults.mat'])
load('Ntrial_type.mat')
load ('evt.mat','evt_trial', 'evt_lick_L', 'evt_lick_R')
load ('time.mat');

% average over trials (100 frames)
% Get_DLC_VAR_Trial_Mat_script

%% OPTIONAL: plot time continuous (raw) DLC-VAR and average (flat)
% plot_continous_DLCvar_script

%% Plot DLC results using accurate trials type selection (corL v corR) or (opto vs omi) ...
% plot_save_real_trialtype_average_GoCueCenter

%% TO compare TTL Lick Sensor vs DLC Lick detection
Plot_DLCLick_vs_TTLSensor_align


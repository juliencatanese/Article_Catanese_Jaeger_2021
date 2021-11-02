% Get_evt_Restrict_TaskEpoch_JCscript
% Create evt restricted to the TASK 
% written by Julien Catanese 2/22/2019
% last updated JC 2/22/2019

clearvars -EXCEPT FolderID MouseID nf

load('evt.mat');
load('time.mat');
load('Epochs_pre_post_task_st_end.mat');

Restrict_StartTask = min(idx_task_st_end);
Restrict_StopTask =max(idx_task_st_end);


time2 = time(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 

evt_trial2 = evt_trial(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 
evt_opto2 = evt_opto(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 
evt_delay2 = evt_delay(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 

evt_lick_L2 = evt_lick_L(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 
evt_lick_R2 = evt_lick_R(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 

evt_puff_R2 = evt_puff_R(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 
evt_puff_L2 = evt_puff_L(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 

evt_rwd_L2 = evt_rwd_L(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH 
evt_rwd_R2 = evt_rwd_R(Restrict_StartTask: Restrict_StopTask); % RESTRICT TO TASK EPOCH

time = time2;

evt_trial = evt_trial2;
evt_opto = evt_opto2;
evt_delay = evt_delay2;

evt_lick_L = evt_lick_L2;
evt_lick_R = evt_lick_R2;

evt_puff_R = evt_puff_R2;
evt_puff_L = evt_puff_L2;

evt_rwd_L = evt_rwd_L2;
evt_rwd_R = evt_rwd_R2;

save('time_Restrict_TaskEpoch.mat', 'time'); 
disp('time saved')

save('evt_Restrict_TaskEpoch.mat', 'evt_trial', 'evt_opto', 'evt_delay',...
    'evt_lick_L', 'evt_lick_R', 'evt_puff_R', 'evt_puff_L', 'evt_rwd_L', 'evt_rwd_R'); 
disp('evt saved')

clearvars -EXCEPT FolderID MouseID nf


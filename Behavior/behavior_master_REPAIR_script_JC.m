%  behavior_master_script_JC
%  play  this for each folder contening "board-DIN-00.dat"
%  "board-DIN-01.dat" ...
%  JCatanese 2017 in  JaegerLab

%%
clear all, close all, 
%% To run a loop over folder 
MouseID = 'JCVGAT17';
FolderID = dir(['D:\JC_Data\Acute\Awake\' MouseID '\w*d*Lv*18*']);

for nf=7:max(size(FolderID))
    cd(['D:\JC_Data\Acute\Awake\' MouseID '\' FolderID(nf).name])
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name])
 
    %% convert .dat into .mat
    if isempty(dir('*.mat'))
    dat2mat_JC_Script
    else
        disp('already computed')
    end
  %%
load('info.mat');
load('evt.mat');
load('time.mat');
sr=20000; 
figure, hold on
plot(time(300*sr:600*sr), evt_trial(300*sr:600*sr),'b')
plot(time(300*sr:600*sr), evt_puff_L(300*sr:600*sr)*0.5,'r')
plot(time(300*sr:600*sr), evt_delay(300*sr:600*sr).*0.3,'g')
legend('trigtrial','puffL','delay')
xlabel('5min time (s)')
saveas(gcf,'evtfig','png')
    
    %% plot the behavioral raster trial/trials with licks    
%     plot_task1_LickRASTER_JC_Script
    
    %% 
%     plot_task1_Perf_Puff_LR_training_JC_Script
    
    %% %% optional not finished
%     Get_task_param_JC_Script
%     plot_task1_Perc_ErrDelay_JC_Script
    
end





%% not finished = need to take delay into account
% Get_Table_trials_behav_JC_Script

%% not finished = need to take delay into account
% plot_task1_ErrResp_Perf_sorted_JC_Script

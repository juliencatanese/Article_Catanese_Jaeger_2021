% Manual_Epochs_St_End_JCspecScript
% Specific Script located in each session folder with different values.
% by JC 11/13/2018
close all;
load info;
MouseID = info.info_notes.MouseID;
Day = info.info_notes.Day;



load('time.mat');
load('evt.mat','evt_opto', 'evt_trial', 'evt_lick_L', 'evt_lick_R', 'evt_puff_L');
evt_lick = evt_lick_L + evt_lick_R;
clear evt_lick_R evt_lick_L;




%% plot All Session evts
figure(1),close(1), h1=figure(1),
hold on,  plot(time, evt_lick*3, 'k');
hold on,  plot(time, evt_trial*2,'r');
hold on,  plot(time, evt_opto, 'c');
hold on,  plot(time, evt_puff_L, 'b');

legend('lick', 'trigTrial', 'opto');
title([MouseID ' ' Day] );

%% MANUALLY DEFINE EPOCHS HERE [idx_st, idx_end] :
disp(' ATTENTION NEED TO CREATE Manual_Epochs_St_End_JCspecScript ' )
disp(' ATTENTION Use the following template to generate it' )
disp(' ATTENTION Replace XXX by time value that you obs in figure1')
disp(' ATTENTION SAVE the new Script as a specScript in the current Session Folder')

idx_preStim_st_end = [find(time==50) find(time==150)]
idx_task_st_end = [find(time==215) find(time==1220)]
idx_postStim_st_end = [find(time==1280 ) find(time==1700)]

%% Calculate Opto stim pulse time in sec
StimSize_pre = median(diff(find(diff(find(evt_opto(idx_preStim_st_end(1):idx_preStim_st_end(2))))>1)))
StimTime_pre = roundn((StimSize_pre/20000),-1) % IN SEC

StimSize_task = median(diff(find(diff(find(evt_opto(idx_task_st_end(1):idx_task_st_end(2))))>1)));
StimTime_task = roundn((StimSize_task/20000),-1) % IN SEC

StimSize_post = median(diff(find(diff(find(evt_opto(idx_postStim_st_end(1):idx_postStim_st_end(2))))>1)));
StimTime_post = roundn((StimSize_post/20000),-1) % IN SEC

%% Plot Epochs Figure
figure(2),close(2), h2 = figure(2),
hold on,  plot(time, evt_opto, 'c')
hold on,  plot(time(idx_preStim_st_end), [0.9 0.9] , 'b', 'LineWidth', 5 )
hold on,  plot(time(idx_task_st_end), [0.9 0.9] , 'k', 'LineWidth', 5 )
hold on,  plot(time(idx_postStim_st_end), [0.9 0.9] , 'r', 'LineWidth', 5 )
ylim([0.8 1])
title([MouseID ' ' Day])
legend('opto pulse', ['PRE (pulse=' num2str(StimTime_pre) 's)'], [ 'TASK (pulse=' num2str(StimTime_task) 's)'] , [ 'POST (pulse=' num2str(StimTime_post) 's)'])
if MouseID == 'vgat11'; legend('opto', ['PRE baseline'], [ 'TASK'] , [ 'POST task BLOCKS']); end

%% SAVING
saveas(h1, 'Session_visu_trigtrials_opto_lick_evt.png')
saveas(h2, 'Epochs_pre_post_task_st_end.png')

saveas(h1, ['D:\JC_Figures\behavior\' MouseID '_' Day 'Session_visu_trigtrials_opto_lick_evt.png'])
saveas(h2, ['D:\JC_Figures\behavior\' MouseID '_' Day 'Epochs_pre_post_task_st_end.png'])

save('Epochs_pre_post_task_st_end.mat','idx_preStim_st_end', 'idx_task_st_end', 'idx_postStim_st_end')
unit='sec'
save('StimSize_Epoch.mat','StimTime_pre', 'StimTime_task', 'StimTime_post', 'unit')
disp('saved')
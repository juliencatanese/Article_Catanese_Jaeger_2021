% plot_Abs_Average_3DLCvar_JCScript
% function plot_Abs_Average_DLCvar_JCfun(VAR1, VAR2, VARname_str1, VARname_str2)
% Julien Catanese 9/23/2020


clc, close all, clear all, 
% FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep';
% MouseID =  'vgat14w14d8';
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW'
MouseID =  'vgat17w10d7'

cd(FileLocation);
load([ MouseID '_DLCresults.mat'])
load('Ntrial_type.mat')
load ('evt.mat','evt_trial', 'evt_lick_L', 'evt_lick_R')
load ('time.mat');

% Get_DLC_VAR_Trial_Mat_script
%% average over trials (100 frames)
NoseY_tr = allTab.NoseY(1:100); TongueY_tr = allTab.TongueY(1:100); WhiskerX_tr = allTab.WhiskerX(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    NoseY_tr = [NoseY_tr allTab.NoseY(itr*100:((itr+1)*100)-1)];
    TongueY_tr = [TongueY_tr allTab.TongueY(itr*100:((itr+1)*100)-1)];
    WhiskerX_tr = [WhiskerX_tr allTab.WhiskerX(itr*100:((itr+1)*100)-1)];    
end

%% cor vs opto trials selection
idx1 = [trial.idx_correct_L trial.idx_correct_R]
idx2 = [trial.idx_correct_L_opto trial.idx_correct_R_opto]
idx_str1 = 'cor'
idx_str2 = 'cor-opto'

VAR_DLC = WhiskerX_tr;
VAR_str = 'WhiskerX'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)

VAR_DLC=TongueY_tr;
VAR_str = 'TongueY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)

VAR_DLC=NoseY_tr;
VAR_str = 'NoseY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)

%% imp vs omiopto trials selection
idx1 = sort([trial.idx_errorDelay_PL_CL, trial.idx_errorDelay_PR_CR, trial.idx_errorDelay_PL_CR,...
    trial.idx_errorDelay_PR_CL])
idx2 = sort([trial.idx_errorDelay_PL_CL_opto, trial.idx_errorDelay_PR_CR_opto,...
    trial.idx_errorDelay_PL_CR_opto, trial.idx_errorDelay_PR_CL_opto])
idx_str1 = 'imp'
idx_str2 = 'imp-opto'

VAR_DLC = WhiskerX_tr;
VAR_str = 'WhiskerX'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)

VAR_DLC=TongueY_tr;
VAR_str = 'TongueY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)

VAR_DLC=NoseY_tr;
VAR_str = 'NoseY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)


%% omi vs omiopto trials selection
idx1 = trial.idx_NoLick
idx2 = trial.idx_NoLick_opto
idx_str1 = 'omi'
idx_str2 = 'omi-opto'

VAR_DLC = WhiskerX_tr;
VAR_str = 'WhiskerX'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)

VAR_DLC=TongueY_tr;
VAR_str = 'TongueY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)

VAR_DLC=NoseY_tr;
VAR_str = 'NoseY'
plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)






% %% Plot VAR using idxL and idxR
% VAR1 = TongueY_tr; VAR_str1 = 'TongueY'; 
% VAR2 = NoseY_tr; VAR_str2 = 'NoseY'; 
% VAR3 = WhiskerX_tr; VAR_str3 = 'WhiskerX'; 
% Thr1 = mean(mean(VAR1(10:30,:)));
% Thr2 = mean(mean(VAR2(10:30,:)));
% Thr3 = mean(mean(VAR3(10:30,:)));
% 
% % trtype='omi'
% if trtype=='imp'
% idx_imp = sort([trial.idx_errorDelay_PL_CL,...
%     trial.idx_errorDelay_PR_CR,...
%     trial.idx_errorDelay_PL_CR,...
%     trial.idx_errorDelay_PR_CL])
% elseif trtype=='cor'
%     idx_imp = trial.idx_correct_R
% elseif trtype=='omi'
%      idx_imp = trial.idx_NoLick
% 
% end
% 
% idx_str = 'cor';
% idx = sort([trial.idx_correct_L trial.idx_correct_R])
% ntr=size(idx,2);
% 
% figure, hold on, 
% K=1
% 
% MeanABS1 = mean(abs(VAR1(:,idx)-Thr1),2); 
% MeanABS2 = mean(abs(VAR2(:,idx)-Thr2),2); 
% MeanABS3 = mean(abs(VAR3(:,idx)-Thr3),2); 
% 
% SemABS1 = (std(abs(VAR1(:,idx)-Thr1)',1)/sqrt(ntr))'; 
% SemABS2 = (std(abs(VAR2(:,idx)-Thr2)',1)/sqrt(ntr))'; 
% SemABS3 = (std(abs(VAR3(:,idx)-Thr3)',1)/sqrt(ntr))'; 
% 
% 
% X=[1:1:100]';
% plot(MeanABS1,'g','LineWidth',2);
% plot(MeanABS2,'c','LineWidth',2);
% plot(MeanABS3,'y','LineWidth',2);
% xline(1.5*25,'k--','LineWidth',2);
% plotshaded(X', [MeanABS1-(SemABS1*K) MeanABS1+(SemABS1*K)],'g');
% plotshaded(X', [MeanABS2-(SemABS2*K) MeanABS2+(SemABS2*K)],'c');
% plotshaded(X', [MeanABS3-(SemABS3*K) MeanABS3+(SemABS3*K)],'y');
% 
% legend(VAR_str1 ,VAR_str2, VAR_str3  ,'GoCue');
% title( ['DLC events' MouseID]);
% 
% saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_LR' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'png')
% saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_LR' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'emf')
% 
% save('Average_ABS_DLC.mat','VAR_str1','VAR_str2','VAR_str3',...
%     'MeanABS1','MeanABS2','MeanABS3',...
%     'SemABS1','SemABS2','SemABS3',...
%     'Thr1','Thr2','Thr3');
% 

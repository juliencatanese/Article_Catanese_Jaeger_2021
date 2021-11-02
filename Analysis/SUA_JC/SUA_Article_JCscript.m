%  SUA_Article_JCscript
% by JC 11-26-2018
%%
clear all, close all,

%% Figure2 Ephys Epochs: SUA examples 
% param 
Center='Delay' % center on the delay start epoch
psth_trial_type = 'cor' % all correct trials
col='k' % black
pre=1500 ; %ms
post= 1500; %ms 

% Example Cell Type1 
cd('D:\JC_Analysis\JCVGAT14\vgat14_w14d2_z4670_15cells_100tr_taskopto_G912_VM_180709_211747_vid_optopost_10mW')
Sh=2; Ch=1; CLUST=1; 
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, Center, psth_trial_type, col, pre, post)

% Example Cell Type2 
cd('D:\JC_Analysis\JCVGAT14\vgat14_w14d2_z4670_15cells_100tr_taskopto_G912_VM_180709_211747_vid_optopost_10mW')
Sh=1; Ch=6; CLUST=2;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, Center, psth_trial_type, col, pre, post)

% Example Cell Type3 
cd('D:\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_30cells_250tr_taskoptoStop215_G912_VM_180729_183133_vid_optopost_11mW')
Sh=4; Ch=2; CLUST=1;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, Center, psth_trial_type, col, pre, post)

%% Figure2 Ephys Epochs: Sequences Matrice

% param 
Center='Delay' % center on the delay start epoch
psth_trial_type = 'cor' % all correct trials
col='k' % black
pre=1500 ; %ms
post= 1500; %ms 

% Population Sequences Matrice
Pop_Seq_DistMatrix_ALLsession_script
% Save 
gcf, 
Ncell=round(max(ylim));
mkdir('D:\JC_Figures\Fig_Article');
saveas(gcf, ['D:\JC_Figures\Fig_Article\Seq_DistMat_NSess' num2str(NSess) '_Ncell' num2str(Ncell)], 'png');
saveas(gcf, ['D:\JC_Figures\Fig_Article\Seq_DistMat_NSess' num2str(NSess) '_Ncell' num2str(Ncell)], 'eps');





